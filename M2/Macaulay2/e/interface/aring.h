#ifndef _aring_h_
#  define _aring_h_

#  include "engine-includes.hpp"

// TODO: fix this
#  if defined(__cplusplus)
class Ring;
class RingElement;
#  else
typedef struct Ring Ring;
typedef struct RingElement RingElement;
#  endif

/**
   ARing interface routines
 */

#  if defined(__cplusplus)
extern "C" {
#  endif

const Ring /* or null */ *rawARingZZp(unsigned long p); /* connected */
/* Expects a prime number p in range 2 <= p <= 32749 */

const Ring /* or null */ *rawARingGaloisFieldM2(
    const RingElement *prim); /* connected */
/* same interface as rawGaloisField, but uses different internal code */

const Ring /* or null */ *rawARingGaloisFieldFlintBig(
    const RingElement *prim); /* connected */
/* same interface as rawGaloisField, but uses Flint GF code designed for
   wordsize p, but too big for lookup tables */

const Ring /* or null */ *rawARingGaloisFieldFlintZech(
    const RingElement *prim); /* connected */
/* same interface as rawGaloisField, but uses Flint GF code designed for
   wordsize p, and uses lookup tables */

const Ring /* or null */ *rawARingTower1(const Ring *K, M2_ArrayString names);
/* create a tower ring with the given variable names and base ring */

const Ring /* or null */ *rawARingTower2(const Ring *R1,
                                         M2_ArrayString new_names);

const Ring /* or null */ *rawARingTower3(const Ring *R1,
                                         engine_RawRingElementArray eqns);

/**
   ARing Flint routines
 */

const Ring /* or null */ *rawARingZZFlint(); /* connected */

const Ring /* or null */ *rawARingQQFlint(); /* connected */

const Ring /* or null */ *rawARingZZpFlint(unsigned long p); /* connected */
/* Expects a prime number p in range 2 <= p <= 2^64-1 */

#  if defined(__cplusplus)
}
#  endif

#endif /* _aring_h_ */

// Local Variables:
// indent-tabs-mode: nil
// End:
