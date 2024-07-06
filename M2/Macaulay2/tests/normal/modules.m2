-- equality of modules
S = QQ[x,y,z]/z^2
assert(S^{-1,1} != S^{1,-1})
assert(image id_(S^1) == S^1)
assert(image matrix"x" != S^{-1})
assert(image matrix"x,y" != image matrix"x")
assert(image matrix"x,y" == image matrix"y,x")
assert(coker matrix"x,y" != coker matrix"x")
assert(coker matrix"x,y" == coker matrix"y,x")
assert(subquotient(matrix"x,y", matrix"x2,y") != subquotient(matrix"x", matrix"x2"))
assert(subquotient(matrix"x,y", matrix"x2,y") == subquotient(matrix"x", matrix"x2,y"))
assert(subquotient(matrix"x,y", matrix"x2,y") != subquotient(matrix"z2", matrix"x2,y"))
assert(subquotient(matrix"x,y", matrix"z2") == image matrix"x,y")
assert(subquotient(id_(S^1), matrix"x,y")   == coker matrix"x,y")

