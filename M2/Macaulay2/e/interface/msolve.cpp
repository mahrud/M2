// Direct in-memory interface to msolve's F4 (libneogb).  See interface/msolve.h.

#include <M2/config.h>

#if HAVE_MSOLVE
// neogb/data.h is a C header with no extern "C" guards of its own, so it has to
// be wrapped here -- but it also includes <gmp.h>, which in this configuration
// pulls in <iosfwd> and other C++ templates that cannot appear inside an
// extern "C" block.  Including data.h's system headers first, at normal
// linkage, means their include guards are already set by the time the wrapped
// include is reached, so only msolve's own declarations end up with C linkage.
#  include <stdint.h>
#  include <stdlib.h>
#  include <stdio.h>
#  include <gmp.h>
#  include <string.h>
#  include <limits.h>
#  include <math.h>
#  ifdef _OPENMP
#    include <omp.h>
#  endif

// msolve's headers come next, and their short macros are undefined right
// afterwards: neogb/data.h defines OFFSET, LENGTH, DEG, MULT, COEFFS and
// friends at global scope, which collide with ordinary identifiers in the M2
// headers included below.  Including msolve first and undefining the macros
// keeps the real prototypes (so the calls below are type checked against the
// installed msolve) without letting the macros leak into M2's code.
extern "C" {
#  include <neogb/f4.h>
#  include <neogb/f4sat.h>
#  include <neogb/res.h>
#  include <neogb/basis.h>
#  include <neogb/engine.h>
#  include <neogb/hash.h>
#  include <neogb/io.h>
}
#  undef OFFSET
#  undef LENGTH
#  undef PRELOOP
#  undef COEFFS
#  undef MULT
#  undef BINDEX
#  undef DEG
#  undef UNROLL
#  undef ORDER_COLUMNS
#  undef PARALLEL_HASHING
#  undef SM_OFFSET
#  undef SM_LEN
#  undef SM_PRE
#  undef SM_CFS
#  undef SM_SIDX
#  undef SM_SMON
#endif

#include "interface/msolve.h"

#include "error.h"
#include "exceptions.hpp"
#include "interface/monomial-ordering.h"
#include "monoid.hpp"
#include "matrices/matrix-con.hpp"
#include "matrices/matrix-stream.hpp"
#include "matrices/matrix.hpp"
#include "monomials/montable.hpp"
#include "rings/geopoly.hpp"
#include "rings/poly.hpp"
#include "rings/polyring.hpp"
#include "rings/ring.hpp"
#include "ring-elements/ring-element.hpp"
#include "resolution/comp-res.hpp"
#include "buffer.hpp"
#include "text-io.hpp"
#include "util.hpp"

#include "../../d/interrupt-jump.h"

#include <algorithm>
#include <cstdlib>
#include <vector>

// Thread-local: see bin/main.cpp.
extern thread_local JumpCell interrupt_jmp;

M2_bool rawMsolvePresent()
{
#if HAVE_MSOLVE
  return true;
#else
  return false;
#endif
}

#if !HAVE_MSOLVE

const Matrix* rawMsolveGB(const Matrix* /*M*/,
                          int /*elim_block_len*/,
                          int /*degree_limit*/,
                          int /*nr_threads*/,
                          int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

const Matrix* rawMsolveModuleGB(const Matrix* /*M*/,
                                int /*module_order*/,
                                int /*degree_limit*/,
                                int /*nr_threads*/,
                                int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

const Matrix* rawMsolveSyzygy(const Matrix* /*M*/,
                              int /*syz_limit*/,
                              int /*syz_rows*/,
                              int /*degree_limit*/,
                              int /*nr_threads*/,
                              int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

const Matrix* rawMsolveSaturate(const Matrix* /*M*/,
                                const Matrix* /*F*/,
                                int /*nr_threads*/,
                                int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

M2_arrayint rawMsolveMinimalBetti(const Matrix* /*M*/,
                                  int /*length_limit*/,
                                  int /*nr_threads*/,
                                  int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

const RingElement* rawMsolvePoincare(const Matrix* /*M*/,
                                     int /*nr_threads*/,
                                     int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

Computation* rawMsolveResolution(const Matrix* /*M*/,
                                 int /*length_limit*/,
                                 int /*nr_threads*/,
                                 int /*info_level*/)
{
  ERROR("this Macaulay2 was built without the msolve library");
  return nullptr;
}

#else

namespace {

// msolve implements exactly one monomial order, standard degree reverse
// lexicographic.  A ring whose GRevLex block carries weights w is still usable,
// because substituting x_i -> x_i^(w_i) carries its order over to msolve's
// exactly: a monomial's weighted degree sum(e_i w_i) becomes the image's total
// degree, and scaling coordinates by positive w_i changes neither which
// coordinate revlex breaks a tie on nor the sign of the difference there.  So
// the substitution is order preserving in both directions, and rather than
// performing it on the polynomials -- which is what callers used to do, at
// interpreter speed, and then undo term by term on the (much larger) result --
// it is applied to the exponents while they are marshalled in and out below.
//
// Returns the weight vector of the GRevLex block, or an empty vector if the
// order is not one msolve can be handed at all.
std::vector<int> grevlexWeights(const MonomialOrdering* mo, int nvars)
{
  std::vector<int> none;
  if (mo == nullptr) return none;
  std::vector<int> wts;
  for (unsigned int i = 0; i < mo->len; i++)
    {
      mon_part p = mo->array[i];
      switch (p->type)
        {
          case MO_GREVLEX:
          case MO_GREVLEX2:
          case MO_GREVLEX4:
            if (not wts.empty()) return none;  // at most one grevlex block
            wts.assign(p->nvars, 1);
            break;
          case MO_GREVLEX_WTS:
          case MO_GREVLEX2_WTS:
          case MO_GREVLEX4_WTS:
            if (not wts.empty()) return none;
            wts.assign(p->nvars, 1);
            if (p->wts != nullptr)
              for (int j = 0; j < p->nvars; j++)
                {
                  // grevlex needs strictly positive weights; a zero or negative
                  // one would make the substitution non-invertible
                  if (p->wts[j] < 1) return none;
                  wts[j] = p->wts[j];
                }
            break;
          case MO_POSITION_UP:
          case MO_POSITION_DOWN:
            break;
          default:
            return none;
        }
    }
  if (static_cast<int>(wts.size()) != nvars) return none;
  return wts;
}

// msolve stores exponents, and the per-block degrees it accumulates from them,
// in a uint16_t.  It never checks for overflow, so scaled exponents that do not
// fit have to be rejected here rather than silently wrapping.
const long msolveMaxExponent = 65535;

// msolve's input is three flat arrays: the number of terms of each generator,
// the exponent vectors concatenated (nvars entries per term, in variable
// order), and the coefficients concatenated.  This is exactly how it stores
// polynomials internally, so building it costs one pass over the matrix.
struct MsolveInput
{
  std::vector<int32_t> lens;   // one entry per generator
  std::vector<int32_t> exps;   // nvars entries per term
  std::vector<int32_t> cfs;    // one entry per term
  std::vector<int32_t> comps;  // one entry per term, 1-based row; module only
};

// appends the columns of M, so that several matrices can be concatenated into
// one input, as msolve's saturation expects
void collectInput(const Matrix* M,
                  const PolyRing* P,
                  const std::vector<int>& wts,
                  MsolveInput& in)
{
  const Ring* KK = P->getCoefficientRing();
  const int nvars = P->n_vars();
  const int ncols = M->n_cols();
  const int32_t charac = static_cast<int32_t>(P->characteristic());

  exponents_t exp = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(nvars));

  Matrix::iterator i(M);
  for (int c = 0; c < ncols; c++)
    {
      int32_t nterms = 0;
      i.set(c);
      for (; i.valid(); i.next())
        {
          Nterm* t = i.entry();
          for (Nterm& s : t)
            {
              P->getMonoid()->to_expvector(s.monom, exp);
              // x_i -> x_i^(w_i), applied to the exponent vector in passing
              long deg = 0;
              for (int j = 0; j < nvars; j++)
                {
                  long e = static_cast<long>(exp[j]) * wts[j];
                  deg += e;
                  if (e > msolveMaxExponent or deg > msolveMaxExponent)
                    throw exc::engine_error(
                        "exponent too large for msolve, which stores exponents "
                        "and degrees in 16 bits");
                  in.exps.push_back(static_cast<int32_t>(e));
                }
              std::pair<bool, long> b = KK->coerceToLongInteger(s.coeff);
              if (not b.first)
                throw exc::engine_error("expected word size coefficients");
              int32_t a = static_cast<int32_t>(b.second);
              if (a < 0) a += charac;
              in.cfs.push_back(a);
              nterms++;
            }
        }
      // msolve rejects a generator with no terms, so zero columns are dropped;
      // they contribute nothing to the ideal anyway.
      if (nterms > 0) in.lens.push_back(nterms);
    }
}

// The basis comes back in the same layout as the input, with bcf an array of
// int32_t coefficients since the characteristic is positive.
// MatrixStream appends terms in the order it receives them, so they have to
// arrive sorted descending in the target ring's order.  With no elimination
// block that is automatic: msolve's degree reverse lexicographic order on the
// scaled exponents is exactly the ring's (possibly weighted) grevlex order, as
// argued above.  With an elimination block msolve orders by the block order
// instead, which is a different order on the same monomials, so the terms have
// to be put back in the ring's order before they are handed over.
void sortTermsDescending(const PolyRing* P,
                         int nterms,
                         int nvars,
                         const std::vector<int32_t>& exps,
                         std::vector<int>& order)
{
  const Monoid* M = P->getMonoid();
  const int msize = M->monomial_size();
  std::vector<int> monoms(static_cast<size_t>(nterms) * msize);
  exponents_t exp = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(nvars));
  for (int t = 0; t < nterms; t++)
    {
      for (int j = 0; j < nvars; j++)
        exp[j] = static_cast<int>(exps[static_cast<size_t>(t) * nvars + j]);
      M->from_expvector(exp, monoms.data() + static_cast<size_t>(t) * msize);
    }
  order.resize(nterms);
  for (int t = 0; t < nterms; t++) order[t] = t;
  std::sort(order.begin(), order.end(), [&](int x, int y) {
    return M->compare(monoms.data() + static_cast<size_t>(x) * msize,
                      monoms.data() + static_cast<size_t>(y) * msize) > 0;
  });
}

const Matrix* buildResult(const PolyRing* P,
                          const std::vector<int>& wts,
                          bool needsSorting,
                          int32_t bld,
                          const int32_t* blen,
                          const int32_t* bexp,
                          const int32_t* bcf)
{
  const int nvars = P->n_vars();
  MatrixStream S(P->make_FreeModule(1));
  int64_t ce = 0, cc = 0;
  std::vector<int32_t> exps;
  std::vector<int> order;
  S.idealBegin(static_cast<size_t>(bld));
  for (int32_t k = 0; k < bld; k++)
    {
      const int32_t nterms = blen[k];
      S.appendPolynomialBegin(static_cast<size_t>(nterms));

      // undoes the substitution: every exponent of a basis element is divisible
      // by its weight, since msolve only ever takes least common multiples and
      // quotients of monomials that came from the scaled input, and those stay
      // in the sublattice spanned by the weights
      exps.resize(static_cast<size_t>(nterms) * nvars);
      for (int32_t t = 0; t < nterms; t++)
        for (int j = 0; j < nvars; j++)
          exps[static_cast<size_t>(t) * nvars + j] = bexp[ce++] / wts[j];

      if (needsSorting)
        sortTermsDescending(P, nterms, nvars, exps, order);
      else
        {
          order.resize(nterms);
          for (int32_t t = 0; t < nterms; t++) order[t] = t;
        }

      for (int32_t i = 0; i < nterms; i++)
        {
          const int t = order[i];
          S.appendTermBegin(0);
          for (int j = 0; j < nvars; j++)
            {
              int32_t e = exps[static_cast<size_t>(t) * nvars + j];
              if (e != 0) S.appendExponent(j, e);
            }
          S.appendTermDone(bcf[cc + t]);
        }
      cc += nterms;
      S.appendPolynomialDone();
    }
  S.idealDone();
  return S.value();
}


// Appends the columns of M as elements of the free module M->rows().  msolve
// considers a smaller component index larger, whereas Macaulay2's Position
// => Up considers a larger row index larger, so components are numbered in
// reverse at this boundary and mapped back in buildModuleResult.
void collectModuleInput(const Matrix* M,
                        const PolyRing* P,
                        const std::vector<int>& wts,
                        MsolveInput& in)
{
  const Ring* KK = P->getCoefficientRing();
  const int nvars = P->n_vars();
  const int ncols = M->n_cols();
  const int32_t charac = static_cast<int32_t>(P->characteristic());

  exponents_t exp = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(nvars));

  Matrix::iterator i(M);
  for (int c = 0; c < ncols; c++)
    {
      int32_t nterms = 0;
      i.set(c);
      for (; i.valid(); i.next())
        {
          // msolve numbers components from 1 and reserves 0 for a plain ring
          // monomial, which is what lets one divisibility test cover both a
          // ring divisor and a divisor in the same component
          const int32_t component =
              static_cast<int32_t>(M->n_rows() - i.row());
          Nterm* t = i.entry();
          for (Nterm& s : t)
            {
              P->getMonoid()->to_expvector(s.monom, exp);
              long deg = 0;
              for (int j = 0; j < nvars; j++)
                {
                  long e = static_cast<long>(exp[j]) * wts[j];
                  deg += e;
                  if (e > msolveMaxExponent or deg > msolveMaxExponent)
                    throw exc::engine_error(
                        "exponent too large for msolve, which stores exponents "
                        "and degrees in 16 bits");
                  in.exps.push_back(static_cast<int32_t>(e));
                }
              std::pair<bool, long> b = KK->coerceToLongInteger(s.coeff);
              if (not b.first)
                throw exc::engine_error("expected word size coefficients");
              int32_t a = static_cast<int32_t>(b.second);
              if (a < 0) a += charac;
              in.cfs.push_back(a);
              in.comps.push_back(component);
              nterms++;
            }
        }
      // a zero column contributes nothing to the submodule, and msolve
      // rejects a generator with no terms
      if (nterms > 0) in.lens.push_back(nterms);
    }
}

// The module counterpart of buildResult.  No sorting is needed: MatrixStream
// buckets terms by component and only requires each bucket to arrive in
// descending order of the ring's monomials, and msolve's basis elements are
// already sorted that way under both of the module orders offered here --
// position over term groups by component and sorts within, and term over
// position sorts globally, of which each component is a subsequence.
const Matrix* buildModuleResult(const FreeModule* F,
                                const std::vector<int>& wts,
                                int32_t bld,
                                const int32_t* blen,
                                const int32_t* bexp,
                                const int32_t* bcomp,
                                const int32_t* bcf,
                                const bool reverseComponents)
{
  const PolyRing* P = F->get_ring()->cast_to_PolyRing();
  const int nvars = P->n_vars();
  MatrixStream S(F);
  int64_t ce = 0, cc = 0;
  S.idealBegin(static_cast<size_t>(bld));
  for (int32_t k = 0; k < bld; k++)
    {
      const int32_t nterms = blen[k];
      S.appendPolynomialBegin(static_cast<size_t>(nterms));
      for (int32_t t = 0; t < nterms; t++)
        {
          // a zero basis element is reported by msolve as a single term with
          // component 0; there is no such row, so skip it entirely
          if (bcomp[cc + t] == 0)
            {
              ce += nvars;
              continue;
            }
          S.appendTermBegin(reverseComponents
                                ? F->rank() - bcomp[cc + t]
                                : bcomp[cc + t] - 1);
          for (int j = 0; j < nvars; j++)
            {
              // undoes the substitution x_i -> x_i^(w_i), as in buildResult
              int32_t e = bexp[ce++] / wts[j];
              if (e != 0) S.appendExponent(j, e);
            }
          S.appendTermDone(bcf[cc + t]);
        }
      cc += nterms;
      S.appendPolynomialDone();
    }
  S.idealDone();
  return S.value();
}

// Everything msolve requires of the ring.  The module entry point shares all
// of it; only the one-row requirement is specific to the ideal ones.
const PolyRing* checkedRing(const Matrix* M)
{
  const PolyRing* P = M->get_ring()->cast_to_PolyRing();
  if (P == nullptr)
    throw exc::engine_error("expected a matrix over a polynomial ring");

  // msolve works over prime fields of word size characteristic only.
  long charac = static_cast<long>(P->characteristic());
  if (charac <= 0 or charac >= (1L << 31))
    throw exc::engine_error(
        "expected a polynomial ring whose characteristic is positive and "
        "less than 2^31");

  // A tower such as (ZZ/p[x,y])[u,v] is a PolyRing whose coefficients are
  // themselves polynomials, which msolve has no way to represent.
  if (P->getCoefficientRing()->cast_to_PolyRing() != nullptr)
    throw exc::engine_error(
        "expected a polynomial ring whose coefficient ring is a field");

  // A quotient S/J is a PolyRing here too, but msolve knows nothing of J and
  // would silently compute in S.  Callers wanting a Groebner basis over a
  // quotient must lift to S and append the presentation of the quotient to the
  // generators; msolveGBMatrix in the Msolve package does exactly that.
  if (P->is_quotient_ring())
    throw exc::engine_error(
        "expected a polynomial ring, not a quotient of one");

  return P;
}

const PolyRing* checkedPolyRing(const Matrix* M)
{
  const PolyRing* P = checkedRing(M);
  if (M->n_rows() != 1)
    throw exc::engine_error("expected a matrix with one row");
  return P;
}

// the grevlex weights of M's ring, refusing any order msolve cannot represent
std::vector<int> checkedWeights(const PolyRing* P)
{
  std::vector<int> wts =
      grevlexWeights(P->getMonoid()->getMonomialOrdering(), P->n_vars());
  if (wts.empty())
    throw exc::engine_error(
        "expected a ring with a degree reverse lexicographic order, possibly "
        "weighted; msolve implements no other monomial order");
  return wts;
}

// The order strategy the resolution entry points run in.
//
// msolve carries this as a res_strat_t with three axes -- the base order on
// R^r, which component index counts as larger, and how the levels above zero
// are derived -- precisely so that it can be varied and measured.  Changing
// what Macaulay2 asks for is changing this one function.
//
// Position over term is the engine's own default and what this has always
// used.  Term over position produces a different (and on modules measurably
// smaller) Groebner basis, frame and differential; only the minimal Betti
// numbers and the Hilbert numerator are invariants across the two, which is
// what msolve's own selftest pins down.  The component direction is left at
// RES_POS_DOWN because collectModuleInput already numbers components in
// reverse to obtain Macaulay2's Position => Up, so asking for UP here as well
// would compose back to Down.
res_strat_t msolveResolutionStrategy()
{
  res_strat_t s = res_strat_default();
  s.base = RES_MORD_POT;
  s.pos = RES_POS_DOWN;
  s.lift = RES_LIFT_SCHREYER;
  return s;
}

// msolve's own defaults, as used by its command line driver.
const int32_t mon_order = 0;  // 0 is degree reverse lexicographic
const int32_t ht_size = 17;
const int32_t max_nr_pairs = 0;
const int32_t reset_ht = 0;
const int32_t la_option = 2;
const int32_t reduce_gb = 1;
const int32_t pbm_file = 0;
const int32_t use_signatures = 0;

void clampOptions(int& nr_threads, int& info_level)
{
  if (nr_threads < 1) nr_threads = 1;
  if (info_level < 0) info_level = 0;
  if (info_level > 2) info_level = 2;
}

// What a Betti table needs of the grading, beyond what a Groebner basis needs
// of the order.
//
// The substitution x_i -> x_i^(w_i) that carries a weighted grevlex ring into
// msolve's unweighted one is harmless to a Groebner basis, because it is order
// preserving and commutes with the lcms and quotients Buchberger takes, so the
// whole computation mirrors itself step for step.  The same argument carries
// the Schreyer frame and the differential across, since those are built from
// lcms and colon quotients too and every monomial involved stays in the
// sublattice spanned by the weights.  So the *ranks* come back right whatever
// the weights are.
//
// The degrees are the part that does not come for free.  With the substitution
// in place msolve grades by the total degree of the substituted exponents,
// which is sum(e_i w_i); a Betti table and a Hilbert numerator are indexed by
// degree, so they mean what they say only when that is the ring's own degree,
// that is when deg(x_i) == w_i.  This used to be checked, and a ring whose
// order carried weights unrelated to its grading was refused, as was a
// multigraded one outright.
//
// Neither is refused any more, because the grading is now handed over instead
// of inferred.  msolve grades by a finitely generated abelian group: every
// entry point takes a res_grading_t giving the degree of each variable, a
// heft, and any torsion factors, computes the Groebner basis in the heft
// degree reverse lexicographic order that grading induces, and reports minimal
// Betti numbers and a Hilbert numerator bucketed by multidegree.  So the
// resolution entry points below marshal the ring's actual degrees across and
// pass the exponents through untouched.
//
// Two things follow.  The monomial order stops mattering to these entry
// points: msolve works in its own heft order and minimal Betti numbers and the
// Hilbert numerator are invariants of the module, not of the order it was
// computed in.  And a weighted ring no longer approaches msolve's 16-bit
// exponent ceiling through its exponents, since nothing multiplies them by the
// weights any more -- only the degree accumulated from them is weighted, and
// that is checked below.
//
// What is still checked is that the grading is one msolve can schedule by: it
// walks the resolution degree by degree, so it needs a heft that is strictly
// positive on every variable, which is the same thing Macaulay2 needs to speak
// of a resolution at all.
struct MsolveGrading
{
  std::vector<int32_t> degs;  // r x nvars, column major, one degree per column
  std::vector<int32_t> heft;  // r
  std::vector<int> hdegs;     // nvars, heft . deg(x_j), for the 16-bit check
  res_grading_t g = {};
};

void collectGrading(const PolyRing* P, MsolveGrading& out)
{
  const Monoid* M = P->getMonoid();
  const Monoid* D = P->degree_monoid();
  const int nvars = P->n_vars();
  const int r = D->n_vars();

  if (r < 1)
    throw exc::engine_error(
        "expected a graded ring; msolve resolves only graded modules");
  // the degshift msolve reports back travels in a fixed size array
  if (r > RES_MTAB_MAXLEN)
    throw exc::engine_error("expected a grading of smaller degree length");

  // Macaulay2's heft vector is exactly what msolve calls one: the linear form
  // on the degree group that the degree by degree schedule counts up in.  A
  // ring can fail to have one -- degrees {1} and {-1}, say -- and then neither
  // system can resolve over it.
  const std::vector<int>& heft = M->get_heft_vector();
  if (static_cast<int>(heft.size()) != r or
      not M->primary_degrees_of_vars_positive())
    throw exc::engine_error(
        "expected a ring with a heft vector positive on every variable; "
        "msolve schedules a resolution degree by degree and would not "
        "terminate without one");

  out.heft.assign(heft.begin(), heft.end());
  out.hdegs = M->primary_degree_of_vars();

  out.degs.resize(static_cast<size_t>(r) * nvars);
  exponents_t de = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(r));
  for (int j = 0; j < nvars; j++)
    {
      D->to_expvector(M->degree_of_var(j), de);
      for (int k = 0; k < r; k++)
        out.degs[static_cast<size_t>(j) * r + k] = static_cast<int32_t>(de[k]);
    }

  // Macaulay2's degree monoid is free, so the grading group is Z^r with no
  // torsion; msolve takes torsion factors here when there are any.
  out.g.r = r;
  out.g.nt = 0;
  out.g.tord = nullptr;
  out.g.degs = out.degs.data();
  out.g.heft = out.heft.data();
}

// The multidegrees of the rows of F, in the order msolve numbers components
// in, which is Macaulay2's reversed -- see collectModuleInput.  Only their
// differences matter, msolve normalizing them internally and reporting the
// shift it applied.
std::vector<int32_t> collectRowDegrees(const FreeModule* F,
                                       const PolyRing* P,
                                       int r)
{
  const Monoid* D = P->degree_monoid();
  const int nrows = F->rank();
  std::vector<int32_t> rowdegs(static_cast<size_t>(nrows) * r);
  exponents_t de = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(r));
  for (int i = nrows - 1, k = 0; i >= 0; i--, k++)
    {
      D->to_expvector(F->degree(i), de);
      for (int c = 0; c < r; c++)
        rowdegs[static_cast<size_t>(k) * r + c] = static_cast<int32_t>(de[c]);
    }
  return rowdegs;
}

// msolve accumulates the heft degree of a monomial in the same 16 bits it
// stores an exponent in.  Without the substitution the exponents are the
// ring's own and collectModuleInput has already checked them; the heft degree
// is the larger of the two under a weighted grading, and is checked here.
void checkHeftDegrees(const MsolveInput& in,
                      int nvars,
                      const std::vector<int>& hdegs)
{
  const size_t nterms = in.exps.size() / static_cast<size_t>(nvars);
  for (size_t t = 0; t < nterms; t++)
    {
      long d = 0;
      for (int j = 0; j < nvars; j++)
        d += static_cast<long>(in.exps[t * nvars + j]) * hdegs[j];
      if (d > msolveMaxExponent)
        throw exc::engine_error(
            "degree too large for msolve, which stores exponents and degrees "
            "in 16 bits");
    }
}

// The one place the grading is still inferred from the order rather than
// handed over: rawMsolveResolution below, whose result is a live
// ResolutionComputation whose free modules Macaulay2 asks for one at a time.
// Those are graded objects and would need a multidegree each -- msolve reports
// them, res_comp_multidegrees does exactly that -- but building them means
// building Macaulay2 FreeModules over the degree monoid rather than an array
// of ints, so that entry point still takes the substitution route and still
// insists the ring's grading be the one the substitution induces.  The one
// shot Betti and Hilbert entry points do not; see collectGrading.
void checkSinglyGradedByWeights(const PolyRing* P, const std::vector<int>& wts)
{
  const Monoid* D = P->degree_monoid();
  if (D->n_vars() != 1)
    throw exc::engine_error(
        "expected a singly graded ring; msolve returns a resolution of a "
        "singly graded module only, minimalBetti aside");

  exponents_t de = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(1));
  for (int j = 0; j < P->n_vars(); j++)
    {
      D->to_expvector(P->getMonoid()->degree_of_var(j), de);
      if (de[0] != wts[j])
        throw exc::engine_error(
            "expected the degree of each variable to equal its weight in the "
            "degree reverse lexicographic block; msolve grades by the "
            "weighted degree, and would index the free modules by it");
    }
}

// What export_module_betti hands back, in both of the indexings it offers.
//
// The heft indexed tables run from degree zero, msolve having normalized the
// row degrees to start there, so degshift is what puts them back on the
// caller's scale.  The multigraded ones are indexed by bucket -- the
// multidegrees that actually occur, which are a sparse subset of a lattice --
// and mdegs holds those degrees with the shift already added back, so they are
// the caller's own.
struct MsolveBetti
{
  int32_t nlevels = 0;
  int32_t maxdeg = 0;
  int32_t degshift = 0;
  std::vector<int32_t> betti;  // nlevels * (maxdeg+1), level major
  std::vector<int32_t> hilb;   // maxdeg+1
  bool trivial = false;        // the submodule is zero, so nothing was run

  int32_t dlen = 1;            // int32 slots per multidegree
  int32_t ndegs = 0;           // number of distinct multidegrees
  std::vector<int32_t> mdegs;  // ndegs * dlen, on the caller's scale
  std::vector<int32_t> mheft;  // ndegs
  std::vector<int32_t> mbetti; // nlevels * ndegs, level major
  std::vector<int32_t> mhilb;  // ndegs
};

// Runs the resolution far enough to fill in what was asked for: with
// wantBetti false only the frame is built and not one field operation happens
// past the Groebner basis, since the Hilbert numerator is the alternating sum
// of the frame ranks.
//
// Returns false, having filled in nothing, if the computation was interrupted;
// bad input throws instead, so that the two are told apart by the caller -- an
// interrupt is the interpreter's to report, not ours.
bool computeBetti(const Matrix* M,
                  int length_limit,
                  bool wantBetti,
                  bool wantHilb,
                  int nr_threads,
                  int info_level,
                  MsolveBetti& out)
{
  const PolyRing* P = checkedRing(M);
  // The order is nothing to these entry points -- msolve computes in the one
  // its grading induces, and what comes back is an invariant of the module --
  // but a ring msolve cannot be handed exponents from at all, a Laurent or
  // local one, is still refused here, and this is where that is decided.
  checkedWeights(P);

  MsolveGrading grading;
  collectGrading(P, grading);

  const int nvars = P->n_vars();
  const int r = grading.g.r;
  const FreeModule* F = M->rows();
  const int nrows = F->rank();
  const long charac = static_cast<long>(P->characteristic());

  if (nrows < 1)
    throw exc::engine_error("expected a matrix with at least one row");
  if (length_limit < 0)
    throw exc::engine_error("expected a non-negative length limit");
  if (wantHilb and length_limit != 0)
    throw exc::engine_error(
        "the Hilbert numerator is an alternating sum over the whole "
        "resolution and admits no length limit");

  // Components are reversed at this boundary, as in collectModuleInput, so
  // the row degrees have to be reversed to match.  Only their differences
  // matter to msolve, which reports the shift it applied.  These are the
  // rows' full multidegrees, r entries each, on the same scale as the degrees
  // of the variables handed over in the grading.
  const std::vector<int32_t> rowdegs = collectRowDegrees(F, P, r);
  const std::vector<int> unitWeights(nvars, 1);

  MsolveInput in;
  collectModuleInput(M, P, unitWeights, in);
  checkHeftDegrees(in, nvars, grading.hdegs);

  // The submodule is zero, so the cokernel is F itself and its minimal
  // resolution is F concentrated in level zero.  msolve rejects an empty
  // input outright, and there is nothing for it to do here anyway.  Both
  // tables are filled in by hand, the heft indexed one from the rows' heft
  // degrees and the multigraded one from their multidegrees.
  if (in.lens.empty())
    {
      std::vector<int32_t> rowheft(nrows);
      for (int i = 0; i < nrows; i++)
        rowheft[i] = static_cast<int32_t>(F->primary_degree(nrows - 1 - i));

      int32_t lo = rowheft[0], hi = rowheft[0];
      for (int i = 1; i < nrows; i++)
        {
          if (rowheft[i] < lo) lo = rowheft[i];
          if (rowheft[i] > hi) hi = rowheft[i];
        }
      out.trivial = true;
      out.nlevels = 1;
      out.degshift = lo;
      out.maxdeg = hi - lo;
      out.betti.assign(static_cast<size_t>(out.maxdeg) + 1, 0);
      out.hilb.assign(static_cast<size_t>(out.maxdeg) + 1, 0);
      out.dlen = r;
      for (int i = 0; i < nrows; i++)
        {
          out.betti[rowheft[i] - lo]++;
          out.hilb[rowheft[i] - lo]++;

          // one bucket per distinct multidegree, found by scanning: a free
          // module has few enough rows for that not to matter
          const int32_t* d = rowdegs.data() + static_cast<size_t>(i) * r;
          int u = 0;
          for (; u < out.ndegs; u++)
            if (std::equal(d, d + r,
                           out.mdegs.begin() + static_cast<size_t>(u) * r))
              break;
          if (u == out.ndegs)
            {
              out.mdegs.insert(out.mdegs.end(), d, d + r);
              out.mheft.push_back(rowheft[i]);
              out.mbetti.push_back(0);
              out.mhilb.push_back(0);
              out.ndegs++;
            }
          out.mbetti[u]++;
          out.mhilb[u]++;
        }
      return true;
    }

  clampOptions(nr_threads, info_level);

  int32_t nlevels = 0;
  int32_t maxdeg = 0;
  int32_t degshift = 0;
  int32_t* betti = nullptr;
  int32_t* hilb = nullptr;
  res_mtable_t mtab = {};

  // msolve does not poll Macaulay2's interrupted flag.  Install the jump
  // target the SIGINT handler uses while control is inside the library, as
  // rawMsolveGB does, so a Ctrl-C during a long resolution comes straight back
  // here.  The jump lands in this frame, so returning from it destroys the
  // input vectors above normally.
  interrupt_jmp.is_set = true;
  if (SETJMP(interrupt_jmp.addr))
    {
      interrupt_jmp.is_set = false;
      return false;
    }

  const res_strat_t resStrat = msolveResolutionStrategy();
  const int64_t nelts = export_module_betti(
      malloc, &nlevels, &maxdeg, &degshift,
      wantBetti ? &betti : nullptr,
      wantHilb ? &hilb : nullptr,
      nullptr, nullptr, nullptr, nullptr,
      // The same Betti numbers and the same Hilbert numerator indexed by the
      // multidegrees that occur rather than collapsed onto the heft degree.
      // Under a singly graded ring the two agree up to the sparseness of the
      // indexing; under a multigraded one this is the only place the
      // multidegrees are to be had, which is what packPoincare needs.
      &mtab,
      // Unlike for a Groebner basis or a syzygy matrix, the strategy is not
      // a user visible choice here: minimal Betti numbers and the Hilbert
      // numerator are invariants of the module, so every strategy gives the
      // same table -- which is exactly what msolve's selftest checks across
      // the whole strategy matrix.  It still changes how long getting there
      // takes, hence msolveResolutionStrategy above.
      in.lens.data(), in.exps.data(), in.comps.data(), in.cfs.data(),
      rowdegs.data(), static_cast<uint32_t>(charac), mon_order, &resStrat,
      // The ring's own grading, so msolve works in the heft degree reverse
      // lexicographic order it induces and indexes what it reports by the
      // ring's degrees.  The exponents above are passed through untouched:
      // it is the grading, not a substitution, that carries the weights now.
      &grading.g,
      static_cast<int32_t>(nvars), static_cast<int32_t>(nrows),
      static_cast<int32_t>(in.lens.size()),
      static_cast<int32_t>(length_limit),
      wantBetti ? 1 : 0,
      0 /* the cheap structural check always runs; the exact d o d = 0 one
         * costs several times the resolution itself */,
      ht_size, static_cast<int32_t>(nr_threads), max_nr_pairs, la_option,
      static_cast<int32_t>(info_level));
  interrupt_jmp.is_set = false;

  // export_module_betti reports every rejection on stderr and returns 0
  // without allocating anything
  if (nelts == 0 or (wantBetti and betti == nullptr)
      or (wantHilb and hilb == nullptr) or mtab.dlen != r
      or (wantHilb and mtab.ndegs > 0 and mtab.hilbnum == nullptr))
    {
      free_module_betti_result_data(free, &betti, &hilb);
      free_module_mtable_data(free, &mtab);
      throw exc::engine_error("msolve returned no Betti table");
    }

  out.nlevels = nlevels;
  out.maxdeg = maxdeg;
  out.degshift = degshift;
  if (wantBetti)
    out.betti.assign(betti,
                     betti + static_cast<size_t>(nlevels) * (maxdeg + 1));
  if (wantHilb)
    out.hilb.assign(hilb, hilb + static_cast<size_t>(maxdeg) + 1);

  // The multigraded tables come back on msolve's normalized scale, degrees
  // having been shifted so the lightest row sits at zero; degshift is added
  // back here so that what leaves this function is in the caller's degrees,
  // as the heft indexed tables are once packBetti adds its own shift.
  out.dlen = mtab.dlen;
  out.ndegs = mtab.ndegs;
  if (mtab.ndegs > 0)
    {
      out.mdegs.assign(
          mtab.degs, mtab.degs + static_cast<size_t>(mtab.ndegs) * mtab.dlen);
      for (int32_t u = 0; u < mtab.ndegs; u++)
        for (int32_t k = 0; k < mtab.dlen; k++)
          out.mdegs[static_cast<size_t>(u) * mtab.dlen + k] +=
              mtab.degshift[k];
      // and likewise for the heft degrees the buckets sit over, degshift
      // being by definition the heft degree of the multidegree just added
      // back -- this is the same shift packBetti applies to the heft indexed
      // table, so that summing the multigraded table over a heft fibre still
      // gives that table's column
      out.mheft.assign(mtab.heft, mtab.heft + mtab.ndegs);
      for (int32_t u = 0; u < mtab.ndegs; u++) out.mheft[u] += degshift;
      if (wantBetti)
        out.mbetti.assign(
            mtab.betti,
            mtab.betti + static_cast<size_t>(mtab.nlevels) * mtab.ndegs);
      if (wantHilb)
        out.mhilb.assign(mtab.hilbnum, mtab.hilbnum + mtab.ndegs);
    }

  free_module_betti_result_data(free, &betti, &hilb);
  free_module_mtable_data(free, &mtab);

  return true;
}

// Packs a Betti table the way unpackEngineBetti reads it: three header
// entries and then one row per slanted degree, degree minus level.  The
// bounds are the extent of the nonzero entries, which is what betti_make
// arrives at for the engine's own resolutions.
M2_arrayint packBetti(const MsolveBetti& b)
{
  const int nlev = b.trivial ? 1 : b.nlevels;
  const int nd = b.maxdeg + 1;

  bool any = false;
  int lo = 0, hi = 0, len = 0;
  for (int l = 0; l < nlev; l++)
    for (int d = 0; d < nd; d++)
      {
        if (b.betti[static_cast<size_t>(l) * nd + d] == 0) continue;
        const int slanted = d + b.degshift - l;
        if (not any)
          {
            lo = hi = slanted;
            any = true;
          }
        else
          {
            if (slanted < lo) lo = slanted;
            if (slanted > hi) hi = slanted;
          }
        if (l > len) len = l;
      }

  const int nrows = hi - lo + 1;
  M2_arrayint result = M2_makearrayint(3 + nrows * (len + 1));
  result->array[0] = lo;
  result->array[1] = hi;
  result->array[2] = len;
  for (int i = 3; i < static_cast<int>(result->len); i++) result->array[i] = 0;
  for (int l = 0; l <= len; l++)
    for (int d = 0; d < nd; d++)
      {
        const int v = b.betti[static_cast<size_t>(l) * nd + d];
        if (v == 0) continue;
        result->array[3 + (d + b.degshift - l - lo) * (len + 1) + l] = v;
      }
  return result;
}

// Packs the multigraded table the way unpackMsolveBetti reads it: the degree
// length and the homological length, then a count of the nonzero entries at
// each level, then those entries themselves, grouped by level, each one a
// multidegree followed by its heft degree and the Betti number.
//
//   [r, len, n_0, ..., n_(len-1), (d_1, ..., d_r, heft, beta) x n_i per level i]
//
// A slanted layout of the kind packBetti uses has nothing to slant by here --
// a multidegree minus a level is not a multidegree -- and the table is far too
// sparse for a dense one: over a grading of length seven the buckets that carry
// anything are a vanishing fraction of the lattice, and each multidegree that
// occurs at all was seen to carry a nonzero entry at exactly one level.  So
// what is sent is the nonzero entries and nothing else, in the level major
// order they are already stored in, which makes the packing a straight scan.
M2_arrayint packMsolveBetti(const MsolveBetti& b)
{
  const int nlev = b.trivial ? 1 : b.nlevels;
  const int r = b.dlen;
  const int stride = r + 2;

  // trailing levels with nothing in them are dropped, as betti_make drops them
  std::vector<int32_t> counts(nlev, 0);
  int len = 0;
  for (int l = 0; l < nlev; l++)
    {
      for (int32_t u = 0; u < b.ndegs; u++)
        if (b.mbetti[static_cast<size_t>(l) * b.ndegs + u] != 0) counts[l]++;
      if (counts[l] != 0) len = l + 1;
    }

  size_t nnz = 0;
  for (int l = 0; l < len; l++) nnz += counts[l];

  M2_arrayint result = M2_makearrayint(2 + len + nnz * stride);
  result->array[0] = r;
  result->array[1] = len;
  for (int l = 0; l < len; l++) result->array[2 + l] = counts[l];

  int* rec = result->array + 2 + len;
  for (int l = 0; l < len; l++)
    for (int32_t u = 0; u < b.ndegs; u++)
      {
        const int32_t v = b.mbetti[static_cast<size_t>(l) * b.ndegs + u];
        if (v == 0) continue;
        for (int k = 0; k < r; k++)
          rec[k] = b.mdegs[static_cast<size_t>(u) * r + k];
        rec[r] = b.mheft[u];
        rec[r + 1] = v;
        rec += stride;
      }
  return result;
}

// The numerator as an element of the degrees ring, ZZ[T_0..T_(r-1)] with the
// T_k inverted, so a generator in a negative degree is fine.
//
// This reads the multigraded table rather than the heft indexed one even when
// r is 1.  The two carry the same information there, but only the multigraded
// one carries it as degrees of the ring rather than as heft degrees, and for
// r > 1 the heft degree does not determine the monomial at all.
const RingElement* packPoincare(const PolyRing* P, const MsolveBetti& b)
{
  const PolynomialRing* D = P->get_degree_ring();
  if (D == nullptr)
    throw exc::engine_error("the ring has no degrees ring");
  const Monoid* DM = D->getMonoid();
  const Ring* ZZring = D->getCoefficientRing();

  if (DM->n_vars() != b.dlen)
    throw exc::engine_error("the degrees ring does not match the grading");

  const PolyRing* DP = D->cast_to_PolyRing();
  if (DP == nullptr)
    throw exc::engine_error("the degrees ring is not a flat polynomial ring");

  // The buckets come out of msolve in whatever order the resolution discovered
  // them in, and every one of them is a distinct monomial here, so the
  // numerator is assembled rather than computed.  Assembling it by summing into
  // one running total is what a term at a time looks like, but Ring::add_to is
  // f = add(f, g) and PolyRing::add copies both of its arguments before merging
  // them, so that spelling copies the whole accumulated polynomial once per
  // term: quadratically many term allocations, which on a module with ten
  // thousand multidegrees cost several times the resolution being reported on
  // and left the run dominated by garbage collection.  A geobucket merges each
  // term across log n levels instead, and destructively, so nothing is copied.
  NtermHeap H(DP);
  monomial mon = DM->make_one();
  exponents_t e = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(b.dlen));
  for (int32_t u = 0; u < b.ndegs; u++)
    {
      if (b.mhilb[u] == 0) continue;
      for (int32_t k = 0; k < b.dlen; k++)
        e[k] = b.mdegs[static_cast<size_t>(u) * b.dlen + k];
      DM->from_expvector(e, mon);
      H.add(D->make_flat_term(ZZring->from_long(b.mhilb[u]), mon));
    }
  ring_elem f = H.value();
  return RingElement::make_raw(D, f);
}

// ---------------------------------------------------------------------------
// A resolution kept alive
//
// Everything else in this file is a single call: marshal in, run msolve,
// marshal out, free.  This one is not, because the question it answers is not
// a single question.  Macaulay2 asks for the shape of a resolution first --
// rawResolutionGetFree wants the rank and the degrees of F_i and nothing else
// -- and only later, and only for the levels it turns out to care about, asks
// for a differential.  Materializing the whole complex to answer the first
// question would defeat the point on exactly the inputs this exists for.
//
// So the object holds an msolve res_comp_t across calls.  Constructing it runs
// the module Groebner basis and the whole Schreyer frame, which is
// combinatorial: after that every free module is known and not one field
// operation has happened past the Groebner basis.  get_matrix is what makes
// msolve reduce, and only up to the level asked for.
//
// The lifetime is the one thing to get right, msolve being malloc/free
// throughout and Macaulay2 being garbage collected.  ResolutionComputation
// derives from our_gc_cleanup, so the destructor -- and with it res_comp_free
// -- runs from a GC finalizer, which intern_res installs below.  Nothing
// msolve allocated is ever handed to the collector: the arrays a differential
// comes back in are copied into Macaulay2 objects and released inside
// get_matrix.
class MsolveResComputation : public ResolutionComputation
{
 public:
  MsolveResComputation(const Matrix* M,
                       const PolyRing* P,
                       std::vector<int> wts,
                       res_comp_t* handle)
      : mInput(M),
        mRing(P),
        mWeights(std::move(wts)),
        mHandle(handle),
        mNLevels(res_comp_nlevels(handle)),
        mDegShift(res_comp_degshift(handle))
  {
    mFrees.assign(mNLevels > 0 ? static_cast<size_t>(mNLevels) : 0, nullptr);
  }

  ~MsolveResComputation() override { res_comp_free(&mHandle); }

  // Whether the frame ended on its own rather than at the length limit.
  bool isComplete() const { return res_comp_is_complete(mHandle) != 0; }

 protected:
  // The only stop condition msolve takes is the length limit, and that was
  // spent when the frame was built.
  bool stop_conditions_ok() override { return true; }

  // The frame is built by the constructor and the differential is built by
  // get_matrix, so there is no third phase for this to drive.
  void start_computation() override {}

  int complete_thru_degree() const override
  {
    throw exc::engine_error(
        "complete_thru_degree is not implemented for msolve resolutions");
  }

  const Matrix* get_matrix(int level) override;
  const FreeModule* get_free(int level) override;
  M2_arrayint get_betti(int type) const override;

  void text_out(buffer& o) const override
  {
    o << "msolve resolution computation, " << mNLevels << " levels";
    if (not isComplete()) o << " (truncated)";
    o << newline;
  }

 private:
  const Matrix* mInput;
  const PolyRing* mRing;
  std::vector<int> mWeights;
  res_comp_t* mHandle;
  int mNLevels;
  int mDegShift;

  // memoized, so that the source of D_i and the target of D_{i+1} are the
  // same object and the complex composes
  std::vector<const FreeModule*> mFrees;
};

const FreeModule* MsolveResComputation::get_free(int level)
{
  if (level < 0 or level >= mNLevels) return mRing->make_FreeModule(0);
  // F_0 is the caller's own ambient free module, in the caller's own row
  // order.  msolve numbers components in reverse (see collectModuleInput), so
  // this is the one level where the two disagree, and get_matrix maps D_1's
  // components back accordingly.
  if (level == 0) return mInput->rows();
  if (mFrees[level] != nullptr) return mFrees[level];

  const int rk = res_comp_rank(mHandle, level);
  if (rk < 0) throw exc::engine_error("msolve lost a level of the resolution");

  std::vector<int32_t> degs(rk > 0 ? static_cast<size_t>(rk) : 1, 0);
  if (rk > 0 and res_comp_degrees(mHandle, level, degs.data()) != 0)
    throw exc::engine_error(
        "msolve could not report the degrees of a free module");

  const Monoid* D = mRing->degree_monoid();
  FreeModule* F = mRing->make_FreeModule();
  monomial deg = D->make_one();
  exponents_t e = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(1));
  for (int i = 0; i < rk; i++)
    {
      // msolve normalized the row degrees to start at zero; degshift is what
      // puts them back on the caller's scale
      e[0] = degs[i] + mDegShift;
      D->from_expvector(e, deg);
      F->append(deg);
    }
  mFrees[level] = F;
  return F;
}

const Matrix* MsolveResComputation::get_matrix(int level)
{
  const FreeModule* tar = get_free(level - 1);
  const FreeModule* src = get_free(level);

  const int rk = (level >= 1 and level < mNLevels) ? src->rank() : 0;
  if (rk == 0)
    {
      MatrixConstructor zero(tar, src);
      return zero.to_matrix();
    }

  int32_t* dlen = nullptr;
  int32_t* dexp = nullptr;
  int32_t* dcomp = nullptr;
  void* dcf = nullptr;

  // msolve does not poll Macaulay2's interrupted flag, so install the jump
  // target the SIGINT handler uses, as everywhere else here.  Unlike the one
  // shot entry points, an interrupt here is recoverable: the levels below
  // rd->thru are complete and untouched, and asking again recomputes from
  // there over the top of whatever the abandoned pass left behind.
  interrupt_jmp.is_set = true;
  if (SETJMP(interrupt_jmp.addr))
    {
      interrupt_jmp.is_set = false;
      return nullptr;
    }
  const int64_t nterms = res_comp_differential(
      malloc, mHandle, level, &dlen, &dexp, &dcomp, &dcf);
  interrupt_jmp.is_set = false;

  if (nterms == 0 or dlen == nullptr or dexp == nullptr or dcomp == nullptr
      or dcf == nullptr)
    {
      free_module_differential_data(free, &dlen, &dexp, &dcomp, &dcf);
      throw exc::engine_error("msolve returned no differential");
    }

  // components of D_1 index F_0, which is the caller's, hence reversed;
  // higher ones index a frame level, which is ours and is not
  const Matrix* raw =
      buildModuleResult(tar, mWeights, rk, dlen, dexp, dcomp,
                        static_cast<const int32_t*>(dcf), level == 1);
  free_module_differential_data(free, &dlen, &dexp, &dcomp, &dcf);

  // The source MatrixStream inferred from the column degrees agrees with
  // get_free(level) entry for entry, but the complex only composes in
  // Macaulay2 if it is literally the free module D_{level+1} targets.
  MatrixConstructor mat(tar, src);
  for (int j = 0; j < rk; j++) mat.set_column(j, raw->elem(j));
  return mat.to_matrix();
}

M2_arrayint MsolveResComputation::get_betti(int type) const
{
  if (type == 4 or type == 0)
    throw exc::engine_error(
        "this resolution is the nonminimal one; ask for minimal Betti numbers "
        "with rawMsolveMinimalBetti, which extracts them from ranks alone and "
        "never materializes a differential");
  if (type != 1)
    throw exc::engine_error(
        "that Betti display is not available for msolve resolutions");

  MsolveBetti b;
  b.nlevels = mNLevels;
  b.degshift = mDegShift;

  std::vector<std::vector<int32_t>> degs(
      mNLevels > 0 ? static_cast<size_t>(mNLevels) : 0);
  for (int l = 0; l < mNLevels; l++)
    {
      const int rk = res_comp_rank(mHandle, l);
      degs[l].assign(rk > 0 ? static_cast<size_t>(rk) : 1, 0);
      if (rk > 0 and res_comp_degrees(mHandle, l, degs[l].data()) != 0)
        throw exc::engine_error("msolve could not report the frame degrees");
      for (int i = 0; i < rk; i++)
        if (degs[l][i] > b.maxdeg) b.maxdeg = degs[l][i];
    }

  const size_t nd = static_cast<size_t>(b.maxdeg) + 1;
  b.betti.assign(static_cast<size_t>(mNLevels) * nd, 0);
  for (int l = 0; l < mNLevels; l++)
    {
      const int rk = res_comp_rank(mHandle, l);
      for (int i = 0; i < rk; i++)
        b.betti[static_cast<size_t>(l) * nd + degs[l][i]]++;
    }
  return packBetti(b);
}

}  // namespace

const Matrix* rawMsolveGB(const Matrix* M,
                          int elim_block_len,
                          int degree_limit,
                          int nr_threads,
                          int info_level)
{
  // rawMsolveGB is the entry point exposed to the interpreter.  Keep the
  // ideal implementation below for one-row matrices, and route genuine
  // module presentations through msolve's module F4 implementation.  The
  // public four-argument interface has no separate module-order option, so
  // use Macaulay2's default: term over position up.
  //
  // TODO: expose a module-order option supporting all four useful variants:
  // position up/down over term and term over position up/down.  msolve
  // currently has POT/TOP with one component direction; reversing component
  // numbering here supplies Up, and the unreversed numbering supplies Down.
  if (M != nullptr and M->n_rows() != 1)
    {
      if (elim_block_len != 0)
        {
          ERROR("msolve does not combine elimination blocks with module "
                "orders");
          return nullptr;
        }
      return rawMsolveModuleGB(M, 1, degree_limit, nr_threads, info_level);
    }

  // A degree ceiling lives on msolve's module entry point and not on its
  // ideal one, and a rank one module basis is the ideal basis element for
  // element -- msolve's own selftest asserts as much -- so a limited request
  // for an ideal goes the module way round.  The module orders coincide
  // there too: with one component there is nothing for them to order.
  if (degree_limit > 0)
    {
      if (elim_block_len != 0)
        {
          ERROR("msolve does not combine elimination blocks with a degree "
                "limit");
          return nullptr;
        }
      return rawMsolveModuleGB(M, 1, degree_limit, nr_threads, info_level);
    }

  try
    {
      const PolyRing* P = checkedPolyRing(M);
      const std::vector<int> wts = checkedWeights(P);
      long charac = static_cast<long>(P->characteristic());
      const int nvars = P->n_vars();
      if (elim_block_len < 0 or elim_block_len > nvars)
        throw exc::engine_error("expected the elimination block length to be "
                                "between 0 and the number of variables");

      MsolveInput in;
      collectInput(M, P, wts, in);

      // With no nonzero generators there is nothing for msolve to do, and it
      // would reject the input outright (calling exit(1) on our behalf).
      if (in.lens.empty())
        {
          MatrixConstructor mat(P->make_FreeModule(1), 0);
          return mat.to_matrix();
        }

      clampOptions(nr_threads, info_level);

      int32_t bld = 0;
      int32_t* blen = nullptr;
      int32_t* bexp = nullptr;
      void* bcf = nullptr;

      // Note msolve reduces in.cfs into [0, charac) in place, which is fine
      // since we own it and do not use it afterwards.
      // msolve does not poll M2's interrupted flag.  Install the jump target
      // used by the SIGINT handler while control is inside the library so a
      // Ctrl-C can return immediately to the interpreter.
      interrupt_jmp.is_set = true;
      if (SETJMP(interrupt_jmp.addr))
        {
          interrupt_jmp.is_set = false;
          return nullptr;
        }
      export_f4(malloc,
                &bld, &blen, &bexp, &bcf,
                in.lens.data(),
                in.exps.data(),
                in.cfs.data(),
                static_cast<uint32_t>(charac),
                mon_order,
                static_cast<int32_t>(elim_block_len),
                static_cast<int32_t>(nvars),
                static_cast<int32_t>(in.lens.size()),
                ht_size,
                static_cast<int32_t>(nr_threads),
                max_nr_pairs,
                reset_ht,
                la_option,
                reduce_gb,
                pbm_file,
                static_cast<int32_t>(info_level));
      interrupt_jmp.is_set = false;

      if (blen == nullptr or bexp == nullptr or bcf == nullptr)
        throw exc::engine_error("msolve returned no basis");

      const Matrix* result =
          buildResult(P, wts, elim_block_len > 0, bld, blen, bexp,
                      static_cast<const int32_t*>(bcf));

      free_f4_julia_result_data(
          free, &blen, &bexp, &bcf, static_cast<int64_t>(bld), charac);

      return result;
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}


const Matrix* rawMsolveModuleGB(const Matrix* M,
                                int module_order,
                                int degree_limit,
                                int nr_threads,
                                int info_level)
{
  try
    {
      const PolyRing* P = checkedRing(M);
      const std::vector<int> wts = checkedWeights(P);
      long charac = static_cast<long>(P->characteristic());
      const int nvars = P->n_vars();
      const FreeModule* F = M->rows();
      const int nrows = F->rank();

      if (nrows < 1)
        throw exc::engine_error("expected a matrix with at least one row");

      // res.h numbers the module orders SCHREYER, POT, TOP.  Components are
      // reversed at the interface boundary, so both choices below have
      // Macaulay2's Position => Up direction.
      //
      // volatile only to quiet -Wclobbered: this is set before the setjmp
      // below and read after it, and the compiler cannot see that the jump
      // path returns without reading it.
      res_strat_t gbStrat = res_strat_default();
      switch (module_order)
        {
          case 0: gbStrat.base = RES_MORD_POT; break;
          case 1: gbStrat.base = RES_MORD_TOP; break;
          default:
            throw exc::engine_error(
                "expected the module order to be 0 (position up over term) "
                "or 1 (term over position up)");
        }

      MsolveInput in;
      collectModuleInput(M, P, wts, in);

      // nothing to do, and msolve would reject an empty input outright
      if (in.lens.empty())
        {
          MatrixConstructor mat(F, 0);
          return mat.to_matrix();
        }

      clampOptions(nr_threads, info_level);

      // No row degrees are handed over, and that is a correctness matter
      // rather than a missed optimization.  row_degs is not just a hint to
      // msolve's degree by degree schedule: set_module_exponent_vector in
      // neogb's res_module.c folds the component's shift into ev[DEG], the
      // very slot cmp_blocks compares first, so the order msolve computes in
      // weighs deg(m) + row_degs[i] before breaking ties in the ring.
      // Macaulay2's order on a free module weighs deg(m) and then the
      // component; the twists of the target play no part in it, and gb of a
      // matrix into S^{0,-1} has the same lead terms as gb of the same
      // entries into S^2.  Passing the twists therefore returns a basis that
      // is a Groebner basis in a different -- perfectly good, but other --
      // order, which msolveDeclareGB's forceGB would declare to be one in
      // Macaulay2's.  Constant row degrees msolve normalizes away, so they
      // were only ever consequential in exactly the case they got wrong.
      //
      // rawMsolveSyzygy and rawMsolveMinimalBetti do hand over a grading, and
      // rightly so: there the twists are part of the answer being asked for,
      // not of an order that has to agree with one Macaulay2 already fixed.
      //
      // The cost is that a twisted target now looks inhomogeneous to msolve
      // and its schedule proceeds by monomial degree alone -- which is what
      // Macaulay2's own engine does on this order in any case.

      // A degree ceiling is only meaningful once msolve and Macaulay2 agree
      // on what a degree is, which is the same condition the row degrees are
      // passed under: with a weighted grevlex block the exponents are
      // substituted before msolve sees them, so its degrees are not the
      // ring's and a ceiling stated in the ring's degrees would cut in the
      // wrong place.  Refusing is the only honest answer -- a silently
      // rescaled ceiling would return the wrong basis without saying so.
      res_stop_t stop = res_stop_none();
      const int32_t deglimit = static_cast<int32_t>(degree_limit);
      if (degree_limit > 0)
        {
          if (not unweighted)
            throw exc::engine_error(
                "msolve cannot honour a degree limit over a ring whose "
                "grevlex block is weighted");
          stop.max_degree = &deglimit;
        }

      int32_t bld = 0;
      int32_t* blen = nullptr;
      int32_t* bexp = nullptr;
      int32_t* bcomp = nullptr;
      void* bcf = nullptr;

      // as in rawMsolveGB, msolve reduces in.cfs in place, which is fine
      // since we own it and do not use it afterwards, and as there it does not
      // poll M2's interrupted flag, so the SIGINT handler's jump target has to
      // be installed around the call
      interrupt_jmp.is_set = true;
      if (SETJMP(interrupt_jmp.addr))
        {
          interrupt_jmp.is_set = false;
          return nullptr;
        }
      int64_t nterms = export_module_f4(malloc,
                &bld, &blen, &bexp, &bcomp, &bcf,
                in.lens.data(),
                in.exps.data(),
                in.comps.data(),
                in.cfs.data(),
                nullptr /* the row degrees; see above */,
                static_cast<uint32_t>(charac),
                mon_order,
                &gbStrat,
                nullptr /* the standard grading; a Groebner basis needs only
                         * the order, and the exponent substitution already
                         * carries a weighted one over -- see collectGrading
                         * for what a graded answer needs beyond that */,
                degree_limit > 0 ? &stop : nullptr,
                static_cast<int32_t>(nvars),
                static_cast<int32_t>(nrows),
                static_cast<int32_t>(in.lens.size()),
                ht_size,
                static_cast<int32_t>(nr_threads),
                max_nr_pairs,
                la_option,
                reduce_gb,
                static_cast<int32_t>(info_level));
      interrupt_jmp.is_set = false;

      // export_module_f4 reports every rejection on stderr and returns 0
      // without allocating anything, rather than calling exit as export_f4
      // does on bad input
      if (nterms == 0 or blen == nullptr or bexp == nullptr
          or bcomp == nullptr or bcf == nullptr)
        throw exc::engine_error("msolve returned no module basis");

      const Matrix* result = buildModuleResult(
          F, wts, bld, blen, bexp, bcomp, static_cast<const int32_t*>(bcf),
          true);

      free_module_f4_result_data(free, &blen, &bexp, &bcomp, &bcf);

      return result;
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const Matrix* rawMsolveSyzygy(const Matrix* M,
                              int syz_limit,
                              int syz_rows,
                              int degree_limit,
                              int nr_threads,
                              int info_level)
{
  try
    {
      const PolyRing* P = checkedRing(M);
      const std::vector<int> wts = checkedWeights(P);
      const int nvars = P->n_vars();
      const int nrows = M->rows()->rank();
      const int ncols = M->cols()->rank();
      const long charac = static_cast<long>(P->characteristic());

      if (nrows < 1)
        throw exc::engine_error("expected a matrix with at least one row");

      // NOTE: unlike rawMsolveGB, which routes modules through Macaulay2's
      // default term over position up order and reproduces gens gb exactly,
      // this returns a *position over term* Groebner basis, and so neither the
      // same elements nor the same column order as gens gb syz.
      //
      // That is not a sorting artifact and cannot be fixed by reordering.  The
      // graph module trick needs an order eliminating the R^nr_rows block, so
      // that the elements whose original part vanishes are exactly the ones
      // whose lead term sits in the adjoined block.  Position over term gives
      // that; term over position does not.  For a homogeneous element of the
      // graph module every term carries the same ev[DEG], so degree reverse
      // lexicographic ties on the degree and picks the lead by reverse lex on
      // the exponents with the component only a tie break -- the lead can land
      // in the adjoined block while the original part is still nonzero, and
      // the characterization collapses.
      //
      // Matching gens gb syz therefore needs msolve to combine an elimination
      // block with a module order, which res.h currently refuses outright.
      MsolveInput in;
      collectModuleInput(M, P, wts, in);
      if (in.lens.empty())
        {
          MatrixConstructor mat(M->cols(), 0);
          return mat.to_matrix();
        }
      if (static_cast<int>(in.lens.size()) != ncols)
        throw exc::engine_error(
            "msolve syzygies do not yet support zero input columns");

      volatile bool unweighted = true;  // volatile as for mord above
      for (int j = 0; j < nvars; j++)
        if (wts[j] != 1) unweighted = false;

      std::vector<int32_t> rowdegs;
      if (unweighted)
        {
          rowdegs.reserve(nrows);
          for (int i = nrows - 1; i >= 0; i--)
            rowdegs.push_back(
                static_cast<int32_t>(M->rows()->primary_degree(i)));
        }

      clampOptions(nr_threads, info_level);

      // The syzygy matrix has one row per column of M, and the components
      // come back unreversed on this path -- see buildModuleResult -- so
      // msolve's first syz_rows components are Macaulay2's first syz_rows
      // rows, with no translation needed.  The degree limit is refused over
      // a weighted grevlex block for the reason rawMsolveModuleGB gives.
      res_stop_t stop = res_stop_none();
      const int32_t deglimit = static_cast<int32_t>(degree_limit);
      if (degree_limit > 0)
        {
          if (not unweighted)
            throw exc::engine_error(
                "msolve cannot honour a degree limit over a ring whose "
                "grevlex block is weighted");
          stop.max_degree = &deglimit;
        }
      if (syz_limit > 0) stop.syz_limit = static_cast<int32_t>(syz_limit);
      if (syz_rows > 0) stop.syz_rows = static_cast<int32_t>(syz_rows);
      const bool limited =
          degree_limit > 0 or syz_limit > 0 or syz_rows > 0;

      int32_t nlevels = 0;
      int32_t* ranks = nullptr;
      int32_t* degs = nullptr;
      int32_t* dlen = nullptr;
      int32_t* dexp = nullptr;
      int32_t* dcomp = nullptr;
      void* dcf = nullptr;

      const res_strat_t resStrat = msolveResolutionStrategy();

      // the resolution does not poll M2's interrupted flag either
      interrupt_jmp.is_set = true;
      if (SETJMP(interrupt_jmp.addr))
        {
          interrupt_jmp.is_set = false;
          return nullptr;
        }
      const int64_t nterms = export_module_resolution(
          malloc, &nlevels, &ranks, &degs, &dlen, &dexp, &dcomp, &dcf,
          in.lens.data(), in.exps.data(), in.comps.data(), in.cfs.data(),
          unweighted ? rowdegs.data() : nullptr,
          static_cast<uint32_t>(charac), mon_order, &resStrat,
          nullptr /* the standard grading, as for rawMsolveModuleGB */,
          limited ? &stop : nullptr,
          static_cast<int32_t>(nvars), static_cast<int32_t>(nrows),
          static_cast<int32_t>(in.lens.size()), 2, RES_SYZ_OF_INPUT,
          0 /* structural verification is intrinsic to the graph module */,
          ht_size, static_cast<int32_t>(nr_threads), max_nr_pairs, la_option,
          static_cast<int32_t>(info_level));
      interrupt_jmp.is_set = false;

      if (nterms == 0 or ranks == nullptr or degs == nullptr
          or dlen == nullptr or dexp == nullptr or dcomp == nullptr
          or dcf == nullptr)
        throw exc::engine_error("msolve returned no syzygy computation");

      const Matrix* result;
      if (nlevels < 3 or ranks[2] == 0)
        {
          MatrixConstructor mat(M->cols(), 0);
          result = mat.to_matrix();
        }
      else
        {
          int64_t firstTerms = 0;
          for (int i = 0; i < ranks[1]; i++) firstTerms += dlen[i];
          result = buildModuleResult(
              M->cols(), wts, ranks[2], dlen + ranks[1],
              dexp + firstTerms * nvars, dcomp + firstTerms,
              static_cast<const int32_t*>(dcf) + firstTerms, false);
        }

      free_module_resolution_result_data(
          free, &ranks, &degs, &dlen, &dexp, &dcomp, &dcf);
      return result;
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const Matrix* rawMsolveSaturate(const Matrix* M,
                                const Matrix* F,
                                int nr_threads,
                                int info_level)
{
  // msolve has no export_f4sat wrapper the way it has export_f4, so this
  // reproduces what its command line driver does for `-S`: build the basis from
  // the ideal generators, build a second basis holding the polynomials to
  // saturate by, and hand both to core_f4sat.
  bs_t* bs = nullptr;
  bs_t* sat = nullptr;
  ht_t* bht = nullptr;
  md_t* st = nullptr;
  try
    {
      const PolyRing* P = checkedPolyRing(M);
      const std::vector<int> wts = checkedWeights(P);
      if (F->get_ring() != M->get_ring())
        throw exc::engine_error("expected both matrices over the same ring");
      if (F->n_rows() != 1)
        throw exc::engine_error("expected a matrix with one row");

      long charac = static_cast<long>(P->characteristic());
      // https://github.com/algebraic-solving/msolve/issues/165
      if (charac < (1L << 16))
        throw exc::engine_error(
            "expected characteristic larger than 2^16 for msolve's F4SAT");

      const int nvars = P->n_vars();
      clampOptions(nr_threads, info_level);
      // msolve's F4SAT corrupts the heap when run with more than one thread in
      // process: a loop of saturations is reliably stable at nr_threads = 1 and
      // reliably aborts with "double free or corruption" at nr_threads = 4,
      // with no freeing of our own involved. Its plain F4 (export_f4, used by
      // rawMsolveGB) has no such problem and is left parallel. msolve's own
      // driver never hits this because it computes once and exits.
      nr_threads = 1;

      // msolve expects the ideal generators followed by the polynomials to
      // saturate by, with nr_nf recording how many of the trailing ones those
      // are; collectInput appends, so the concatenation is built in place.
      MsolveInput in;
      collectInput(M, P, wts, in);
      const int32_t ngens = static_cast<int32_t>(in.lens.size());
      collectInput(F, P, wts, in);
      const int32_t nsat = static_cast<int32_t>(in.lens.size()) - ngens;

      if (nsat == 0)
        throw exc::engine_error("expected a nonzero polynomial to saturate by");
      if (ngens == 0)
        {
          MatrixConstructor mat(P->make_FreeModule(1), 0);
          return mat.to_matrix();
        }

      int success = initialize_gba_input_data(&bs, &bht, &st,
          in.lens.data(), in.exps.data(), in.cfs.data(),
          static_cast<uint32_t>(charac), mon_order, 0 /* elim_block_len */,
          static_cast<int32_t>(nvars), ngens + nsat, nsat /* nr_nf */,
          ht_size, static_cast<int32_t>(nr_threads), max_nr_pairs, reset_ht,
          la_option, use_signatures, reduce_gb, pbm_file,
          0 /* truncate_lifting */, static_cast<int32_t>(info_level));

      // every generator was zero, so the saturation is the zero ideal
      if (success == -1)
        {
          MatrixConstructor mat(P->make_FreeModule(1), 0);
          return mat.to_matrix();
        }
      if (success == 0) throw exc::engine_error("msolve rejected the input");

      // initialize_gba_input_data forces 32 bit internals via the modulus it is
      // given; the true characteristic is recorded separately, exactly as
      // msolve's driver does.
      st->gfc = static_cast<len_t>(charac);

      sat = initialize_basis(st, bht);
      import_input_data(sat, st, ngens, ngens + nsat,
                        in.lens.data(), in.exps.data(), in.cfs.data(), NULL);
      sat->ld = sat->lml = nsat;
      for (int32_t k = 0; k < nsat; k++) sat->lmps[k] = k;

      int error = 0;
      // Like export_f4, core_f4sat does not poll M2's interrupted flag.
      interrupt_jmp.is_set = true;
      if (SETJMP(interrupt_jmp.addr))
        {
          interrupt_jmp.is_set = false;
          return nullptr;
        }
      success = core_f4sat(bs, sat, st, &error);
      interrupt_jmp.is_set = false;
      // core_f4sat releases sat itself, via free_basis_elements followed by
      // free_basis_without_hash_table, so it must not be freed again below.
      sat = nullptr;
      if (not success or error)
        throw exc::engine_error("msolve failed to saturate");

      int32_t bld = 0;
      int32_t* blen = nullptr;
      int32_t* bexp = nullptr;
      void* bcf = nullptr;
      int64_t nterms = export_results_from_gba(
          &bld, &blen, &bexp, &bcf, malloc, &bs, &bht, &st);
      if (nterms == 0 or blen == nullptr or bexp == nullptr or bcf == nullptr)
        throw exc::engine_error("msolve returned no basis");

      const Matrix* result =
          buildResult(P, wts, false, bld, blen, bexp,
                      static_cast<const int32_t*>(bcf));

      free_f4_julia_result_data(
          free, &blen, &bexp, &bcf, static_cast<int64_t>(bld), charac);

      // same cleanup export_f4 performs after export_results_from_f4
      free_shared_hash_data(bht);
      free_basis(&bs);
      free(st);
      return result;
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawMsolveMinimalBetti(const Matrix* M,
                                  int length_limit,
                                  int nr_threads,
                                  int info_level)
{
  try
    {
      MsolveBetti b;
      if (not computeBetti(M, length_limit, true /* betti */,
                           false /* hilbert */, nr_threads, info_level, b))
        return nullptr;  // interrupted; the interpreter reports it
      return packMsolveBetti(b);
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement* rawMsolvePoincare(const Matrix* M,
                                     int nr_threads,
                                     int info_level)
{
  try
    {
      const PolyRing* P = checkedRing(M);
      MsolveBetti b;
      if (not computeBetti(M, 0 /* no length limit */, false /* betti */,
                           true /* hilbert */, nr_threads, info_level, b))
        return nullptr;  // interrupted; the interpreter reports it
      return packPoincare(P, b);
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

Computation* rawMsolveResolution(const Matrix* M,
                                 int length_limit,
                                 int nr_threads,
                                 int info_level)
{
  try
    {
      const PolyRing* P = checkedRing(M);
      const std::vector<int> wts = checkedWeights(P);
      // the free modules of a resolution are graded objects, so this needs
      // what a graded answer needs of the grading, not merely what a Groebner
      // basis needs of the order
      checkSinglyGradedByWeights(P, wts);

      const int nvars = P->n_vars();
      const int nrows = M->rows()->rank();
      const long charac = static_cast<long>(P->characteristic());

      if (nrows < 1)
        throw exc::engine_error("expected a matrix with at least one row");
      if (length_limit < 0)
        throw exc::engine_error("expected a nonnegative length limit");

      MsolveInput in;
      collectModuleInput(M, P, wts, in);
      if (in.lens.empty())
        throw exc::engine_error(
            "expected a nonzero matrix; msolve does not resolve the zero "
            "submodule");

      std::vector<int32_t> rowdegs;
      rowdegs.reserve(nrows);
      for (int i = nrows - 1; i >= 0; i--)
        rowdegs.push_back(static_cast<int32_t>(M->rows()->primary_degree(i)));

      clampOptions(nr_threads, info_level);

      // volatile only to silence -Wclobbered; C11 makes an object indeterminate
      // after longjmp when it is modified *between* setjmp and longjmp, which
      // this is not
      const res_strat_t resStrat = msolveResolutionStrategy();
      res_comp_t* volatile handle = nullptr;

      interrupt_jmp.is_set = true;
      if (SETJMP(interrupt_jmp.addr))
        {
          interrupt_jmp.is_set = false;
          return nullptr;  // the interpreter reports the interrupt
        }
      handle = res_comp_new(
          in.lens.data(), in.exps.data(), in.comps.data(), in.cfs.data(),
          rowdegs.data(), static_cast<uint32_t>(charac), mon_order,
          // as for the one shot resolution: the Schreyer order the
          // differential runs in is the one position over term induces
          &resStrat,
          nullptr /* the standard grading, as for rawMsolveModuleGB */,
          static_cast<int32_t>(nvars),
          static_cast<int32_t>(nrows),
          static_cast<int32_t>(in.lens.size()),
          static_cast<int32_t>(length_limit), ht_size,
          static_cast<int32_t>(nr_threads), max_nr_pairs, la_option,
          static_cast<int32_t>(info_level));
      interrupt_jmp.is_set = false;

      if (handle == nullptr)
        throw exc::engine_error("msolve could not start the resolution");

      ResolutionComputation* C = new MsolveResComputation(M, P, wts, handle);
      // hands the handle to the collector: remove_res deletes C, and ~C is
      // what calls res_comp_free
      intern_res(C);
      return C;
    } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

#endif  // HAVE_MSOLVE

// Restricts a Groebner basis G of <m> + J*S^r, computed over the ambient ring
// S of a quotient R0 = S/J (by rawMsolveGB/rawMsolveModuleGB on m0 | rels, as
// msolveGBMatrix in the Msolve package does), to a Groebner basis of <m> over
// R0 itself.  Unlike the rest of this file, this needs nothing from msolve: it
// is plain engine bookkeeping, so it is compiled regardless of HAVE_MSOLVE.
//
// A column survives exactly when its lead term -- in the sense of
// leadTerm(Matrix), i.e. the single term of the vector that is largest in the
// free module order, at whatever row it falls in -- is not divisible by any
// lead term of J.  Which row it falls in does not matter: J's initial ideal is
// the same in every row of S^r, since rels is J's presentation tensored with
// the identity.  R0 already has that initial ideal as a MonomialTable, built
// once when R0 was constructed for its own ring arithmetic (normal_form), so
// this is a membership query against existing data, not a second Groebner
// basis computation of leadTerm(rels) the way the Msolve package used to do it
// at the interpreter level.
//
// Note the quotient ideal's MonomialTable stores its generators at component
// 1, not 0 -- see QRingInfo_field::QRingInfo_field in e/rings/qring.cpp.
const Matrix* rawMsolveGBRestrictToQuotient(const Matrix* G, const Ring* R0)
{
  try
    {
      const PolyRing* P = G->get_ring()->cast_to_PolyRing();
      if (P == nullptr)
        throw exc::engine_error("expected a matrix over a polynomial ring");
      const PolynomialRing* PR0 = R0->cast_to_PolynomialRing();
      if (PR0 == nullptr)
        throw exc::engine_error("expected a polynomial ring");

      MonomialTable* ringtable = PR0->get_quotient_MonomialTable();
      const Monoid* Mo = P->getMonoid();
      exponents_t exp = ALLOCATE_EXPONENTS(EXPONENT_BYTE_SIZE(P->n_vars()));

      const Matrix* LT = G->lead_term(-1);
      Matrix::iterator i(LT);
      std::vector<int> keep;
      for (int c = 0; c < G->n_cols(); c++)
        {
          i.set(c);
          if (not i.valid()) continue;  // a zero column is already accounted for
          Nterm* t = i.entry();
          Mo->to_expvector(t->monom, exp);
          if (ringtable == nullptr or ringtable->find_divisor(exp, 1) < 0)
            keep.push_back(c);
        }
      return G->sub_matrix(stdvector_to_M2_arrayint(keep));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// Local Variables:
// indent-tabs-mode: nil
// End:
