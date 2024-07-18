R = QQ[x,y]
p = map(R,QQ)
f = matrix {{x-y, x+2*y, 3*x-y}};
assert( kernel f === image map(R^{{-1},{-1},{-1}},R^{{-1},{-2}},{{-7, -x-2*y}, {-2, x-y}, {3, 0}}) )
g = map(R^1,QQ^3,p,f)
assert not isHomogeneous g
assert( kernel g == image map(QQ^3,QQ^1,{{-7}, {-2}, {3}}) )
assert( coimage g == cokernel map(QQ^3,QQ^1,{{-7}, {-2}, {3}}) )
g = map(R^1,,p,f,Degree => {1})
assert isHomogeneous g
R = QQ[x, Degrees => {{2:1}}]
M = R^1
S = QQ[z];
N = S^1
p = map(R,S,{x},DegreeMap => x -> join(x,x))
assert isHomogeneous p
f = matrix {{x^3}}
g = map(M,N,p,f,Degree => {3,3})
assert isHomogeneous g
assert (kernel g == 0)
assert( coimage g == S^1 )

end
-- FIXME
S = QQ[x,y]
f = basis(1, image 1_S)
target f
map(target f, , f)

F = image f
debug Core

restart
gbTrace=1
errorDepth=1
R = QQ[a..d];
M = coker matrix {{a*b, a*c, b*c, a*d}}
N = image map(R^{4:-2}, , {{-c, -c, 0, -d}, {b, 0, -d, 0}, {0, a, 0, 0}, {0, 0, c, b}})
P = image map(R^{4:-3}, , {{d}, {0}, {-c}, {b}})
f = basis(5, P)
F = image f
inducedMap(P, , gens F) -- 4 gbTrace lines


f = vars R
inducedMap(image f, , f)

f = basis(4, map(P, R^1, {P_0})) -- TODO: this should work
f = basis(4, map(P, R^1, matrix P_0))
B = image f
inducedMap(P, , f) -- good
inducedMap(B, , f) -- bad?

f = basis(5, P)
F = image f
target f
F.Monomials
gbTrace = 3
inducedMap(target f, , F.Monomials)
inducedMap(target f, F)
map(target f, F, F.Monomials
target f == F.Monomials.target
F.Monomials


gens F
cover f
