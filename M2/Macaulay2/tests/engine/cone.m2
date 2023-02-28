debug Core

A = matrix{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}; B = map(ZZ^7,ZZ^0,0);
o = transpose map(ZZ, rawFourierMotzkinEqs(raw transpose A, raw transpose B))

(A', B') = (o_{0..o_(0,numcols o-1)-1}, o_{o_(0,numcols o-1)..numcols o-2})
o = transpose map(ZZ, rawFourierMotzkinEqs(raw transpose A', raw transpose B'))

assert(A == hermite o_{0..o_(0,numcols o-1)-1})
assert(B == o_{o_(0,numcols o-1)..numcols o-2})

needsPackage "FourierMotzkin"
assert((A, B) == fourierMotzkin fourierMotzkin(A, B))
