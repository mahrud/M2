// See BUILD/build/BarvinokFeature.md for the design of this interface.

#include "interface/barvinok.h"

#include <M2/config.h>

#if HAVE_BARVINOK

#  include <M2/math-include.h>

#  include "debug.hpp"
#  include "interface/matrix.h"
#  include "matrix-con.hpp"
#  include "matrix.hpp"
#  include "relem.hpp"

#  include <csetjmp>
#  include <csignal>
#  include <gmpxx.h>
#  include <vector>

// PolyLib's Matrix (a plain typedef, not namespaced) would otherwise clash
// with M2's own engine::Matrix; rename just that one token while pulling in
// barvinok/isl/PolyLib. (Wrapping the includes in a C++ namespace instead
// doesn't work here: barvinok.h's C++-only tail pulls in NTL's QQ, colliding
// with the nested std::less lookup as well as M2's own `using QQ =
// mpq_srcptr` in ringelem.hpp. We only need barvinok_enumerate_with_options,
// declared in that same extern "C" block, so declare it by hand below
// instead of including barvinok/barvinok.h and its C++-only gen_fun tail.)
#  define Matrix BarvinokMatrix
#  include <barvinok/evalue.h>
#  include <barvinok/options.h>
#  include <barvinok/util.h>
#  include <isl/aff.h>
#  include <isl/ctx.h>
#  include <isl/polynomial.h>
#  include <isl/space.h>
#  include <isl/val.h>
#  include <isl/val_gmp.h>
#  undef Matrix

extern "C" {
evalue *barvinok_enumerate_with_options(Polyhedron *P, Polyhedron *C,
                                         struct barvinok_options *options);
}

namespace {

// Build a PolyLib constraint Matrix directly from an M2 Matrix's integer
// entries, in the "combined constraint" layout: column 0 is the equality(0)
// or inequality(1) flag, followed by coefficients, followed by the constant.
// TODO: Check that M is over ZZ (get_mpz() on a non-ZZ ring_elem, e.g. QQ, is
// undefined behavior); callers currently must lift to ZZ themselves.
BarvinokMatrix *toPolyLibMatrix(const Matrix *M)
{
  const size_t r = M->n_rows();
  const size_t c = M->n_cols();
  BarvinokMatrix *mat = Matrix_Alloc(r, c);
  for (size_t i = 0; i < r; i++)
    for (size_t j = 0; j < c; j++)
      mpz_set(mat->p[i][j], M->elem(i, j).get_mpz());
  return mat;
}

// One floor()/div term's affine argument: coeff_0..coeff_{nparam-1}, const,
// all over the single shared integer denom, i.e. floor((coeff.t + const)/denom).
struct DivDef
{
  std::vector<mpz_class> coeff;
  mpz_class cst;
  mpz_class denom;
};

// One monomial of a chamber's quasipolynomial: (coeffNum/coeffDen) * prod_i
// t_i^exp[i] * prod_j floor(divdefs[j])^divExp[j] (divExp[j] == 0 means the
// j-th div is unused in this term).
struct Term
{
  mpz_class coeffNum, coeffDen;
  std::vector<long> exp;
  std::vector<long> divExp;
};

struct ChamberData
{
  const Matrix *facets;
  std::vector<DivDef> divdefs;
  std::vector<Term> terms;
};

// isl_val -> mpz_class, asserting it is an integer (caller already arranged
// this, e.g. by multiplying a rational coefficient by its shared denominator).
mpz_class valToInt(isl_val *v)
{
  mpz_class z;
  isl_val_get_num_gmp(v, z.get_mpz_t());
  isl_val_free(v);
  return z;
}

struct TermCollectorData
{
  unsigned nparam;
  isl_size ndivs;
  ChamberData *chamber;
  bool haveDivDefs;
};

isl_stat collectTerm(__isl_take isl_term *term, void *user)
{
  auto *ud = static_cast<TermCollectorData *>(user);
  Term t;
  t.coeffNum = 0;
  t.coeffDen = 1;
  {
    isl_val *c = isl_term_get_coefficient_val(term);
    isl_val_get_num_gmp(c, t.coeffNum.get_mpz_t());
    isl_val_get_den_gmp(c, t.coeffDen.get_mpz_t());
    isl_val_free(c);
  }
  for (unsigned i = 0; i < ud->nparam; i++)
    t.exp.push_back(isl_term_get_exp(term, isl_dim_param, i));
  for (isl_size j = 0; j < ud->ndivs; j++)
    t.divExp.push_back(isl_term_get_exp(term, isl_dim_div, j));

  // div definitions are shared by every term of the same qpolynomial, so
  // only need to be extracted once (from whichever term arrives first).
  if (!ud->haveDivDefs)
    {
      for (isl_size j = 0; j < ud->ndivs; j++)
        {
          isl_aff *aff = isl_term_get_div(term, j);
          isl_val *denom = isl_aff_get_denominator_val(aff);
          DivDef d;
          d.denom = valToInt(isl_val_copy(denom));
          for (unsigned i = 0; i < ud->nparam; i++)
            {
              isl_val *c = isl_aff_get_coefficient_val(aff, isl_dim_param, i);
              d.coeff.push_back(valToInt(isl_val_mul(c, isl_val_copy(denom))));
            }
          {
            isl_val *cst = isl_aff_get_constant_val(aff);
            d.cst = valToInt(isl_val_mul(cst, isl_val_copy(denom)));
          }
          isl_val_free(denom);
          isl_aff_free(aff);
          ud->chamber->divdefs.push_back(std::move(d));
        }
      ud->haveDivDefs = true;
    }

  ud->chamber->terms.push_back(std::move(t));
  isl_term_free(term);
  return isl_stat_ok;
}

// For some degenerate chambers, isl_qpolynomial_from_evalue hits a real
// barvinok/isl bug: evalue_isl.c's assertion `e->x.p->type == polynomial ||
// flooring || fractional` fails and calls abort() (see the FIXME on the
// disabled benchmark assertion in Chambers.m2 for a concrete repro). An
// assert()-triggered abort() is not a C++ exception -- it can't be caught
// by try/catch, unlike everything else this function might throw -- so the
// only way to keep it from taking down the whole M2 process is to catch the
// resulting SIGABRT with a signal handler and jump back out via siglongjmp.
// This is scoped to SIGABRT specifically, not SIGSEGV: an assert() failure
// happens at a clean, controlled point before any memory is corrupted, so
// recovering from it and continuing is reasonably safe; a genuine SIGSEGV
// would mean actual corruption already happened, and longjmp'ing past that
// would just continue running atop already-broken state.
//
// sigaction's disposition is process-wide, not per-thread, so there's a
// narrow race if some unrelated abort() fires on another thread while this
// one is inside a guarded region: that thread's handler invocation sees
// barvinokRecoveryActive false (correctly, since it's thread_local) and
// falls through to SIG_DFL+raise, which momentarily changes the *global*
// disposition and could in principle race with this thread's guard being
// torn down. Accepted as out of scope: aborting at all is already a rare,
// degenerate-input path, and two unrelated aborts racing is rarer still.
thread_local sigjmp_buf barvinokRecoveryPoint;
thread_local volatile std::sig_atomic_t barvinokRecoveryActive = 0;

extern "C" void barvinokAbortHandler(int sig)
{
  if (barvinokRecoveryActive) siglongjmp(barvinokRecoveryPoint, sig);
  // not inside a guarded region -- restore default behavior and re-raise,
  // same as if this handler had never been installed
  signal(sig, SIG_DFL);
  raise(sig);
}

struct SignalGuard
{
  struct sigaction oldAction{};
  SignalGuard()
  {
    struct sigaction newAction{};
    newAction.sa_handler = barvinokAbortHandler;
    sigemptyset(&newAction.sa_mask);
    sigaction(SIGABRT, &newAction, &oldAction);
    barvinokRecoveryActive = 1;
  }
  ~SignalGuard()
  {
    barvinokRecoveryActive = 0;
    sigaction(SIGABRT, &oldAction, nullptr);
  }
};

// Extract the chamber (a cone in the parameter space) as an M2 Matrix of
// facet inequalities, directly from the Polyhedron's GMP constraint data;
// equality rows are dropped (only ever used the inequality rows of a
// chamber, matching the earlier PolyLib-text-based implementation).
const Matrix *chamberFacets(const Ring *R, const Polyhedron *D)
{
  const unsigned dim = D->Dimension;
  size_t nrows = 0;
  for (unsigned i = 0; i < D->NbConstraints; i++)
    if (mpz_sgn(D->Constraint[i][0]) != 0) nrows++;

  MatrixConstructor mat(R->make_FreeModule(nrows), dim);
  size_t row = 0;
  for (unsigned i = 0; i < D->NbConstraints; i++)
    {
      if (mpz_sgn(D->Constraint[i][0]) == 0) continue;
      for (unsigned j = 0; j < dim; j++)
        mat.set_entry(row, j, R->from_int(D->Constraint[i][j + 1]));
      row++;
    }
  return mat.to_matrix();
}

}  // namespace

const Matrix /* or null */ *rawBarvinokEnumerate(const Matrix *M, const Matrix *C)
{
  try
    {
      SignalGuard signalGuard;
      if (sigsetjmp(barvinokRecoveryPoint, 1) != 0)
        throw exc::engine_error(
            "barvinok/isl aborted while decomposing chambers or converting "
            "a quasipolynomial (a known bug for some degenerate chambers; "
            "see BUILD/build/BarvinokFeature.md)");

      const Ring *R = M->get_ring();

      struct barvinok_options *options = barvinok_options_new_with_defaults();
      isl_ctx *ctx = isl_ctx_alloc_with_options(&barvinok_options_args, options);
      options->chambers = BV_CHAMBERS_ISL;

      BarvinokMatrix *Pm = toPolyLibMatrix(M);
      Polyhedron *P = Constraints2Polyhedron(Pm, options->MaxRays);
      Matrix_Free(Pm);

      BarvinokMatrix *Cm = toPolyLibMatrix(C);
      Polyhedron *Ctx = Constraints2Polyhedron(Cm, options->MaxRays);
      Matrix_Free(Cm);

      evalue *EP = barvinok_enumerate_with_options(P, Ctx, options);
      Enumeration *en = partition2enumeration(EP);

      std::vector<ChamberData> chambers;
      unsigned nparam = 0;
      isl_size maxDivs = 0;
      for (Enumeration *node = en; node != nullptr; node = node->next)
        {
          nparam = node->ValidityDomain->Dimension;
          ChamberData cd;
          cd.facets = chamberFacets(R, node->ValidityDomain);

          isl_space *space = isl_space_params_alloc(ctx, nparam);
          // see the SignalGuard comment above: this can abort() for some
          // degenerate chambers, which is caught around the whole function
          isl_qpolynomial *qp = isl_qpolynomial_from_evalue(space, &node->EP);
          isl_size ndivs = isl_qpolynomial_dim(qp, isl_dim_div);
          if (ndivs > maxDivs) maxDivs = ndivs;

          TermCollectorData ud{nparam, ndivs, &cd, false};
          isl_qpolynomial_foreach_term(qp, collectTerm, &ud);
          isl_qpolynomial_free(qp);

          chambers.push_back(std::move(cd));
        }

      Enumeration_Free(en);  // this also frees the data reachable from EP;
                             // calling evalue_free(EP) afterwards double-frees
                             // (the top-level `evalue` struct itself is not
                             // freed here, and leaks a few bytes per call --
                             // partition2enumeration doesn't expose a way to
                             // release just the outer struct safely)
      Polyhedron_Free(P);
      Polyhedron_Free(Ctx);
      isl_ctx_free(ctx);  // also frees options, which it took ownership of

      // row = [tag, chamberIndex, index, data...]; data is zero-padded on
      // the right to a common width, wide enough for the widest row kind
      // (a term row, 2 + nparam + maxDivs entries -- see barvinok.h).
      const size_t dataWidth = 2 + nparam + maxDivs;
      const size_t totalCols = 3 + dataWidth;

      size_t totalRows = 0;
      for (const auto &cd : chambers)
        totalRows += cd.facets->n_rows() + cd.divdefs.size() + cd.terms.size();

      MatrixConstructor mat(R->make_FreeModule(totalRows), totalCols);
      size_t row = 0;
      for (size_t ci = 0; ci < chambers.size(); ci++)
        {
          const ChamberData &cd = chambers[ci];
          for (size_t i = 0; i < cd.facets->n_rows(); i++)
            {
              mat.set_entry(row, 0, R->from_long(0));  // tag: facet
              mat.set_entry(row, 1, R->from_long(ci));
              mat.set_entry(row, 2, R->from_long(i));
              for (size_t j = 0; j < nparam; j++)
                mat.set_entry(row, 3 + j, cd.facets->elem(i, j));
              row++;
            }
          for (size_t i = 0; i < cd.divdefs.size(); i++)
            {
              const DivDef &d = cd.divdefs[i];
              mat.set_entry(row, 0, R->from_long(1));  // tag: divdef
              mat.set_entry(row, 1, R->from_long(ci));
              mat.set_entry(row, 2, R->from_long(i));
              for (size_t j = 0; j < nparam; j++)
                mat.set_entry(row, 3 + j, R->from_int(d.coeff[j].get_mpz_t()));
              mat.set_entry(row, 3 + nparam, R->from_int(d.cst.get_mpz_t()));
              mat.set_entry(row, 3 + nparam + 1, R->from_int(d.denom.get_mpz_t()));
              row++;
            }
          for (size_t i = 0; i < cd.terms.size(); i++)
            {
              const Term &t = cd.terms[i];
              mat.set_entry(row, 0, R->from_long(2));  // tag: term
              mat.set_entry(row, 1, R->from_long(ci));
              mat.set_entry(row, 2, R->from_long(i));
              mat.set_entry(row, 3, R->from_int(t.coeffNum.get_mpz_t()));
              mat.set_entry(row, 4, R->from_int(t.coeffDen.get_mpz_t()));
              for (size_t j = 0; j < nparam; j++)
                mat.set_entry(row, 5 + j, R->from_long(t.exp[j]));
              for (size_t j = 0; j < t.divExp.size(); j++)
                mat.set_entry(row, 5 + nparam + j, R->from_long(t.divExp[j]));
              row++;
            }
        }
      return mat.to_matrix();
  } catch (const exc::engine_error &e)
    {
      ERROR(e.what());
      return nullptr;
  }
}

#else  // !HAVE_BARVINOK

const Matrix /* or null */ *rawBarvinokEnumerate(const Matrix *M, const Matrix *C)
{
  ERROR("barvinok library not found; reconfigure cmake with barvinok "
        "(and isl, polylib) discoverable via pkg-config to use this function");
  return nullptr;
}

#endif
