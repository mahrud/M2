#ifndef M2_INTERFACE_TORIC_H_
#define M2_INTERFACE_TORIC_H_

#include "engine-includes.hpp"
#include "m2-types.h"

#ifdef __cplusplus
class Matrix;
extern "C" {
#else
typedef struct Matrix Matrix;
#endif

/**
 * Refine a fan given by its rays and maximal cones.
 *
 * The result is a ZZ matrix with tagged rows.  Tag 0 rows are
 * [0, dimension, ray coordinates], and tag 1 rows are
 * [1, number-of-rays, ray indices].  Rows are zero padded to a
 * common width.  All input rays are emitted first and unchanged.
 *
 * `strategy` selects the pulling-order reversal and the Normaliz local mode,
 * `seed` selects the global order perturbation, `limit` is a soft
 * maximal-cone limit (zero means unlimited), and `threads` and `verbose`
 * control the engine computation.
 */
const Matrix *rawSimplicialFan(const Matrix *rays, M2_arrayint cones,
                               int strategy, int seed, int limit,
                               int threads, bool verbose);

/**
 * Refine a fan to a smooth fan.  The strategy selects the pulling-order
 * reversal and Normaliz local mode, with the same remaining option meanings
 * as rawSimplicialFan.  The result uses the same tagged-row encoding.
 */
const Matrix *rawSmoothFan(const Matrix *rays, M2_arrayint cones,
                           int strategy, int seed, int limit,
                           int threads, bool verbose);

#ifdef __cplusplus
}
#endif

#endif
