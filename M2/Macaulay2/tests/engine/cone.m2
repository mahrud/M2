debug Core

(A, B) = (map(ZZ^2, ZZ^4, {{1, -3, 1, 0}, {0, 1, 0, 1}}), map(ZZ^2, ZZ^0, 0))
o = transpose map(ZZ, rawFourierMotzkin(raw transpose A, raw transpose B))

(A, B) = (map(ZZ^3, ZZ^4, 0) || id_(ZZ^4), map(ZZ^7, ZZ^0, 0))
o = transpose map(ZZ, rawFourierMotzkin(raw transpose A, raw transpose B))

(A', B') = (o_{0..o_(0,numcols o-1)-1}, o_{o_(0,numcols o-1)..numcols o-2})
o = transpose map(ZZ, rawFourierMotzkin(raw transpose A', raw transpose B'))

assert(A == hermite o_{0..o_(0,numcols o-1)-1})
assert(B == o_{o_(0,numcols o-1)..numcols o-2})

needsPackage "FourierMotzkin"
assert((A, B) == fourierMotzkin fourierMotzkin(A, B))

(A, B) = (map(ZZ^5, ZZ^3, 0), map(ZZ^5, ZZ^0, 0))
assert((B, id_(ZZ^5)) == fourierMotzkin(A, B))
assert((B, B) == fourierMotzkin fourierMotzkin(A, B))

assert(findHeft {{-1, 0}, {-3, 1}, {-1, 0}, {0, 1}} == {-1, 1})
assert(findHeft {{ 1, 0}, {-3, 1}, { 1, 0}, {0, 1}} == { 1, 4})
