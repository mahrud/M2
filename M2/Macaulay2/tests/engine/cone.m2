debug Core

(A, B) = (map(ZZ^4, ZZ^2, {{1, 0}, {-3, 1}, {1, 0}, {0, 1}}), map(ZZ^4, ZZ^0, 0))

(A, B) = (map(ZZ^3, ZZ^4, 0) || id_(ZZ^4), map(ZZ^7, ZZ^0, 0))
o = transpose map(ZZ, rawFourierMotzkin(raw transpose A, raw transpose B))

(A', B') = (o_{0..o_(0,numcols o-1)-1}, o_{o_(0,numcols o-1)..numcols o-2})
o = transpose map(ZZ, rawFourierMotzkin(raw transpose A', raw transpose B'))

assert(A == hermite o_{0..o_(0,numcols o-1)-1})
assert(B == o_{o_(0,numcols o-1)..numcols o-2})

needsPackage "FourierMotzkin"
assert((A, B) == fourierMotzkin fourierMotzkin(A, B))
