-- from http://nemocas.org/benchmarks.html
R = ZZFlint[x,y,z,t] -- R = ZZ[x,y,z,t]
f = 1 + x + y + z + t
p = f^30;
print "ready"
time q = p*(p + 1);
print "done"
print q

-- not using flint for multiplication, see:
 3# fmpz_set_ui at /home/mahrud/Projects/M2/M2/M2/submodules/flint2/fmpz.h:229
 4# M2::ARingZZ::from_ring_elem(long&, ring_elem const&) const at /home/mahrud/Projects/M2/M2/M2/BUILD/build/../../Macaulay2/e/aring-zz-flint.hpp:232
 5# PolyRing::mult_by_term(ring_elem, ring_elem, int const*) const at /home/mahrud/Projects/M2/M2/M2/BUILD/build/../../Macaulay2/e/poly.cpp:768
 6# PolyRing::mult(ring_elem, ring_elem) const at /home/mahrud/Projects/M2/M2/M2/BUILD/build/../../Macaulay2/e/poly.cpp:802
 7# RingElement::operator*(RingElement const&) const at /home/mahrud/Projects/M2/M2/M2/BUILD/build/../../Macaulay2/e/relem.cpp:84
 8# engine_star__1 at /home/mahrud/Projects/M2/M2/M2/Macaulay2/d/engine.dd:119
 9# actors_star_ at /home/mahrud/Projects/M2/M2/M2/Macaulay2/d/actors.d:241
