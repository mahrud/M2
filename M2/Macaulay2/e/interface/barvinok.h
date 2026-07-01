#ifndef _barvinok_h_
#  define _barvinok_h_

#  include "engine-includes.hpp"

#  if defined(__cplusplus)
class Matrix;
#  else
typedef struct Matrix Matrix;
#  endif

/**
   Vector partition function interface routines (via barvinok/isl)
 */

#  if defined(__cplusplus)
extern "C" {
#  endif

/**************************************************/
/**** Chamber decomposition (via barvinok) ********/
/**************************************************/

// M encodes the data+parameter polytope and C the parameter/context polytope,
// both in the PolyLib "combined constraint" format: one row per constraint,
// first column 0 (equality) or 1 (inequality), followed by the coefficients,
// followed by the constant term. This is the same layout Chambers.m2 used to
// hand-format as text for the barvinok_enumerate CLI's stdin.
//
// Encodes all chambers of the secondary fan, and the (piecewise) quasi-
// polynomial on each, as a single ZZ matrix of tagged rows -- entirely
// numeric, with no string round-trip through isl's printer. Let nparam be
// the number of columns of C minus 2 (i.e. the ambient dimension of the
// chambers), and let D be (numcols of the result) - 5 - nparam (the largest
// number of div/floor terms appearing in any monomial, across all chambers).
// Each row is [tag, chamberIndex, index, data...], zero-padded on the right
// to a common width of 3 + max(nparam, nparam+2, 2+nparam+D):
//   tag 0 (facet):  data = facet coefficients (nparam values); the row is
//                   one inequality of the chamber's cone, in the parameter
//                   space (index is just the row's position, unused).
//   tag 1 (divdef): data = [coeff_0..coeff_{nparam-1}, const, denom] for
//                   the affine argument of the index-th floor()/div term
//                   used by this chamber's quasipolynomial.
//   tag 2 (term):   data = [coeffNum, coeffDen, exp_0..exp_{nparam-1},
//                   divExp_0..divExp_{D-1}]: one monomial of the chamber's
//                   quasipolynomial, i.e. (coeffNum/coeffDen) * prod_i
//                   t_i^exp_i * prod_j floor(divdef_j)^divExp_j (a divExp
//                   of 0 means that div is unused in this term; index is
//                   just the term's position, unused).
const Matrix /* or null */ *rawBarvinokEnumerate(const Matrix *M, const Matrix *C);

#  if defined(__cplusplus)
}
#  endif

#endif /* _barvinok_h_ */
