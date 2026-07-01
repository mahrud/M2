// See BUILD/build/BarvinokFeature.md for the design of this interface.

#include "interface/barvinok.h"

#include <M2/config.h>

#if HAVE_BARVINOK

#  include <M2/math-include.h>

#  include "debug.hpp"
#  include "interface/gmp-util.h"
#  include "interface/matrix.h"
#  include "matrix-con.hpp"
#  include "matrix.hpp"
#  include "relem.hpp"
#  include "util.hpp"

#  include <regex>
#  include <string>
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
#  include <isl/ctx.h>
#  include <isl/polynomial.h>
#  include <isl/printer.h>
#  include <isl/space.h>
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

// Extract the chamber (a cone in the parameter space) as an M2 Matrix of
// facet inequalities, directly from the Polyhedron's GMP constraint data;
// equality rows are dropped, matching the previous PolyLib-text-based
// readPolyLibCone (which only ever used the inequality rows of a chamber).
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
        {
          mpz_ptr z = newitem(__mpz_struct);
          mpz_init_set(z, D->Constraint[i][j + 1]);
          mpz_reallocate_limbs(z);
          mat.set_entry(row, j, ring_elem(z));
        }
      row++;
    }
  return mat.to_matrix();
}

// Wrap every occurrence of floor(...) in an extra pair of parentheses, e.g.
// floor((x)/3)^2 -> (floor((x)/3))^2, matching the M2-side fixup that used
// to be applied by Chambers.m2's readISLQuasipolynomial via a recursive
// PCRE pattern (std::regex has no recursive subpatterns, so this is done
// with a manual balanced-parenthesis scan instead).
std::string parenthesizeFloors(const std::string &s)
{
  std::string out;
  out.reserve(s.size() + 8);
  const std::string floorTok = "floor(";
  size_t i = 0;
  while (i < s.size())
    {
      if (s.compare(i, floorTok.size(), floorTok) == 0)
        {
          size_t j = i + floorTok.size();
          int depth = 1;
          while (j < s.size() && depth > 0)
            {
              if (s[j] == '(') depth++;
              else if (s[j] == ')') depth--;
              j++;
            }
          out += '(';
          out.append(s, i, j - i);
          out += ')';
          i = j;
        }
      else
        {
          out += s[i];
          i++;
        }
    }
  return out;
}

// isl prints most products with an explicit '*' (e.g. "5/2 + x", "3/2 * y^2"),
// but affine expressions nested inside floor() arguments can come out as
// e.g. "2x"; M2 has no implicit multiplication, so insert the '*' back in.
std::string insertImplicitStars(const std::string &s)
{
  static const std::regex re("([0-9])([a-zA-Z])");
  return std::regex_replace(s, re, "$1*$2");
}

// Converts a per-chamber evalue quasipolynomial into an M2 function-literal
// string "(t0,...,tn) -> EXPR", ready to be handed to M2's `value`.
std::string quasipolynomialToM2String(isl_ctx *ctx,
                                       unsigned nparam,
                                       evalue *EP)
{
  isl_space *space = isl_space_params_alloc(ctx, nparam);
  std::vector<std::string> names(nparam);
  for (unsigned i = 0; i < nparam; i++)
    {
      names[i] = "t" + std::to_string(i);
      space = isl_space_set_dim_name(space, isl_dim_param, i, names[i].c_str());
    }

  // NOTE: for some degenerate chambers this hits a real barvinok/isl bug --
  // evalue_isl.c's assert(e->x.p->type == polynomial || flooring || fractional)
  // fails and aborts the whole process (not a catchable C++ exception). See
  // the FIXME on the disabled benchmark assertion in Chambers.m2 for a
  // concrete repro. Likely fix: print the evalue natively instead of
  // routing through isl's qpolynomial conversion.
  isl_qpolynomial *qp = isl_qpolynomial_from_evalue(space, EP);
  isl_printer *pr = isl_printer_to_str(ctx);
  pr = isl_printer_print_qpolynomial(pr, qp);
  char *cstr = isl_printer_get_str(pr);
  std::string printed(cstr);
  free(cstr);
  isl_printer_free(pr);
  isl_qpolynomial_free(qp);

  // printed looks like "[t0, t1] -> { EXPR }"; take the body between the
  // outermost braces so we can prepend our own, already-known parameter list.
  size_t open = printed.find('{');
  size_t close = printed.rfind('}');
  std::string body = (open != std::string::npos && close != std::string::npos && close > open)
                          ? printed.substr(open + 1, close - open - 1)
                          : printed;

  body = parenthesizeFloors(body);
  body = insertImplicitStars(body);

  // trim whitespace
  size_t a = body.find_first_not_of(" \t\n");
  size_t b = body.find_last_not_of(" \t\n");
  body = (a == std::string::npos) ? std::string("0") : body.substr(a, b - a + 1);

  std::string header = "(";
  for (unsigned i = 0; i < nparam; i++)
    {
      if (i > 0) header += ",";
      header += names[i];
    }
  header += ") -> ";
  return header + body;
}

}  // namespace

engine_RawMatrixStringPairArrayOrNull rawBarvinokEnumerate(const Matrix *M,
                                                            const Matrix *C)
{
  try
    {
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

      std::vector<const Matrix *> facetsList;
      std::vector<std::string> polyList;
      for (Enumeration *node = en; node != nullptr; node = node->next)
        {
          const unsigned nparam = node->ValidityDomain->Dimension;
          facetsList.push_back(chamberFacets(R, node->ValidityDomain));
          polyList.push_back(quasipolynomialToM2String(ctx, nparam, &node->EP));
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

      const size_t n = facetsList.size();
      engine_RawMatrixStringPairArray result =
          getmemarraytype(engine_RawMatrixStringPairArray, n);
      result->len = static_cast<int>(n);
      for (size_t i = 0; i < n; i++)
        {
          engine_RawMatrixStringPair pair = new engine_RawMatrixStringPair_struct;
          pair->a = facetsList[i];
          pair->b = string_std_to_M2(polyList[i]);
          result->array[i] = pair;
        }
      return result;
  } catch (const exc::engine_error &e)
    {
      ERROR(e.what());
      return nullptr;
  }
}

#else  // !HAVE_BARVINOK

engine_RawMatrixStringPairArrayOrNull rawBarvinokEnumerate(const Matrix *M,
                                                            const Matrix *C)
{
  ERROR("barvinok library not found; reconfigure cmake with barvinok "
        "(and isl, polylib) discoverable via pkg-config to use this function");
  return nullptr;
}

#endif
