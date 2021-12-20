// Copyright 2014 Michael E. Stillman

#include "aring-gf-flint.hpp"
#include "poly.hpp"
#include "relem.hpp"
#include "ringmap.hpp"
#include "util.hpp"

namespace M2 {

ARingGFFlint::ARingGFFlint(const PolynomialRing& R, const ring_elem a)
    : mOriginalRing(&R),
      mCharacteristic(R.characteristic()),
      mGeneratorComputed(false)
{
  ring_elem b = R.copy(a);
  mPrimitiveElement = RingElement::make_raw(&R, b);
  ring_elem minpoly = mOriginalRing->quotient_element(0);

  std::vector<long> poly;
  RingElement F(&R, minpoly);
  F.getSmallIntegerCoefficients(poly);

  // warning: flint expects unsigned longs, so we must make
  // all of these coeffs non-negative.
  for (long i = poly.size() - 1; i >= 0; i--)
    {
      long a = poly[i];
      if (a == 0) continue;
      if (a < 0) poly[i] += characteristic();
    }

  mDimension = poly.size() - 1;

#if 0
    printf("minpoly: ");
    for (long i=0; i<poly.size(); i++)
      printf(" %ld", poly[i]);
    printf("\n");
#endif
  nmod_poly_t mMinPoly;
  nmod_poly_init(mMinPoly, mCharacteristic);

  for (long i = poly.size() - 1; i >= 0; i--)
    if (poly[i] != 0) nmod_poly_set_coeff_ui(mMinPoly, i, poly[i]);

  fq_nmod_ctx_init_modulus(mBigContext, mMinPoly, "a");
  fq_zech_ctx_init_fq_nmod_ctx(mContext, mBigContext);
#if 0
    fq_zech_ctx_print(mContext);
#endif
  nmod_poly_clear(mMinPoly);

  flint_randinit(mRandomState);
}

ARingGFFlint::ARingGFFlint(const long p0, const long d)
    : mOriginalRing(0),
      mCharacteristic(p0),
      mDimension(d),
      mGeneratorComputed(false)
{
  fmpz_t p;
  fmpz_init_set_ui(p, p0);
  fq_nmod_ctx_init(mBigContext, p, d, "a");
  fq_zech_ctx_init_fq_nmod_ctx(mContext, mBigContext);
#if 0
  fq_zech_ctx_print(mContext);
#endif
  flint_randinit(mRandomState);
}

ARingGFFlint::~ARingGFFlint()
{
  fq_zech_ctx_clear(mContext);
  fq_nmod_ctx_clear(mBigContext);
  mPrimitiveElement = 0;
  flint_randclear(mRandomState);

  if (mGeneratorComputed) fq_zech_clear(&mCachedGenerator, mContext);
}

void ARingGFFlint::getSmallIntegerCoefficients(const ElementType& a,
                                               std::vector<long>& poly) const
{
  fq_nmod_t f;
  fq_nmod_init(f, mBigContext);
  fq_zech_get_fq_nmod(f, &a, mContext);
  long deg = nmod_poly_degree(f);
  poly.resize(deg + 1);
  for (long i = deg; i >= 0; i--) poly[i] = nmod_poly_get_coeff_ui(f, i);
  fq_nmod_clear(f, mBigContext);
}

void ARingGFFlint::fromSmallIntegerCoefficients(
    ElementType& result,
    const std::vector<long>& poly) const
{
  fq_nmod_t f;
  fq_nmod_init(f, mBigContext);

#if 0
    printf("input = ");
    for (long i=0; i<poly.size(); i++)
      printf("%ld ", poly[i]);
    printf("\n");
#endif
  for (long i = poly.size() - 1; i >= 0; i--)
    {
      long a = poly[i];
      if (a == 0) continue;
      if (a < 0) a += characteristic();
      nmod_poly_set_coeff_ui(f, i, a);
    }
#if 0
    printf("  result before reduction = ");
    fq_nmod_print_pretty(f, mBigContext);
    printf("\n");
#endif
  fq_nmod_reduce(f, mBigContext);
#if 0
    printf("  result = ");
    fq_nmod_print_pretty(f, mBigContext);
    printf("\n");
#endif
  fq_zech_set_fq_nmod(&result, f, mContext);
#if 0
    printf("  zech result = %lu", result.value);
    printf("\n");
#endif
  fq_nmod_clear(f, mBigContext);
}

void ARingGFFlint::getGenerator(ElementType& result_gen) const
{
  if (not mGeneratorComputed)
    {
      fq_zech_init(&mCachedGenerator, mContext);
      fq_zech_gen(&mCachedGenerator, mContext);
      mGeneratorComputed = true;
    }
  copy(result_gen, mCachedGenerator);
}

M2_arrayint ARingGFFlint::getModulusPolynomialCoeffs() const
{
  nmod_poly_t modulus;
  nmod_poly_init(modulus, mCharacteristic);
  nmod_poly_set(modulus, fq_zech_ctx_modulus(mContext));
  long deg = nmod_poly_degree(modulus);
  std::vector<long> poly(deg+1);
  for (long i = deg; i >= 0; i--) poly[i] = nmod_poly_get_coeff_ui(modulus, i);
  nmod_poly_clear(modulus);
  return stdvector_to_M2_arrayint(poly);
}

bool ARingGFFlint::promote(const Ring* Rf,
                           const ring_elem f,
                           ElementType& result) const
{
  // ZZ/p[x]/F(x) ---> GF(p,n) = Rf
  if (mOriginalRing != Rf) return false;
  std::vector<long> poly;
  RingElement F(Rf, f);
  F.getSmallIntegerCoefficients(poly);
  fromSmallIntegerCoefficients(result, poly);
  return true;
}

void ARingGFFlint::lift_to_original_ring(ring_elem& result,
                                         const ElementType& f) const
{
  // This code needs review, and tests.  See git issue #612
  std::vector<long> poly;
  getSmallIntegerCoefficients(f, poly);
  result =
      originalRing().getNumeratorRing()->fromSmallIntegerCoefficients(poly, 0);
}

bool ARingGFFlint::lift(const Ring* Rg,
                        const ElementType& f,
                        ring_elem& result) const
{
  // Rg = ZZ/p[x]/F(x) ---> GF(p,n)
  if (mOriginalRing != Rg) return false;
  lift_to_original_ring(result, f);
  return true;
}

void ARingGFFlint::eval(const RingMap* map,
                        const ElementType& f,
                        int first_var,
                        ring_elem& result) const
{
  // f is represented by: f.value, the power of the generator
  if (is_zero(f))
    result = map->get_ring()->zero();
  else
    result =
        map->get_ring()->power(map->elem(first_var), static_cast<int>(f.value));
}

void ARingGFFlint::text_out(buffer& o) const
{
  o << "GF(" << mCharacteristic << "^" << mDimension << ",Flint)";
}

void ARingGFFlint::elem_text_out(buffer& o,
                                 const ElementType& a,
                                 bool p_one,
                                 bool p_plus,
                                 bool p_parens) const
{
  o << fq_zech_get_str_pretty(&a, mContext);
}

int ARingGFFlint::compare_elems(const ElementType& f,
                                const ElementType& g) const
{
  if (f.value > g.value) return GT;
  if (f.value < g.value) return LT;
  return EQ;
}
};

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
