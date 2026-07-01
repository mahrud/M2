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
// Returns one (facets, quasipolynomial) pair per chamber of the secondary fan:
// facets is the ZZ matrix of the chamber's cone inequalities (in the
// parameter space), and quasipolynomial is an M2 function-literal string of
// the form "(t0,...,tn) -> EXPR" ready to be handed to M2's `value`, e.g.
// "(t0,t1) -> ((1 + t0) + (5/2 + t0) * t1 + 3/2 * t1^2)".
engine_RawMatrixStringPairArrayOrNull rawBarvinokEnumerate(const Matrix *M,
                                                            const Matrix *C);

#  if defined(__cplusplus)
}
#  endif

#endif /* _barvinok_h_ */
