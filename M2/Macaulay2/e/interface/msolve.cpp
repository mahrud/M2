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
#include "rings/poly.hpp"
#include "rings/ring.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>

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
  std::vector<int32_t> lens;  // one entry per generator
  std::vector<int32_t> exps;  // nvars entries per term
  std::vector<int32_t> cfs;   // one entry per term
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

// Everything msolve requires of the ring, checked once for all entry points.
const PolyRing* checkedPolyRing(const Matrix* M)
{
  const PolyRing* P = M->get_ring()->cast_to_PolyRing();
  if (P == nullptr)
    throw exc::engine_error("expected a matrix over a polynomial ring");
  if (M->n_rows() != 1)
    throw exc::engine_error("expected a matrix with one row");

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

}  // namespace

const Matrix* rawMsolveGB(const Matrix* M,
                          int elim_block_len,
                          int nr_threads,
                          int info_level)
{
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
      success = core_f4sat(bs, sat, st, &error);
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

#endif  // HAVE_MSOLVE

// Local Variables:
// indent-tabs-mode: nil
// End:
