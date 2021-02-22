Y = Proj (QQ[x])
HH^0 ( OO_Y , Degree => 0 )
G = OO_Y
M = module G
e = 0
M = M / saturate 0_M
A = ring M;
F = presentation A
R = ring F
N = coker lift(presentation M,R) ** coker F
r = numgens R
wR = R^{-r}
E1 = Ext^(r-1)(N,wR)
degrees E1
min degrees E1

--
k = ZZ/101
R = k[a,b,c,d]/(a^4+b^4+c^4+d^4)
X = Proj R
result = table(3, 3, (p, q) -> timing ((p, q) => rank HH^q(cotangentSheaf(p, X))))
assert( {{1, 0, 1}, {0, 20, 0}, {1, 0, 1}} === applyTable(result, last @@ last) )
print new MatrixExpression from result

-- Example 4.1: the bounds can be sharp.
S = QQ[w,x,y,z];
X = Proj S;
I = monomialCurveIdeal(S, {1,3,4})
N = S^1/I;
assert(Ext^1(OO_X, N~(>= 0)) == prune truncate(0, Ext^1(truncate(2, S^1), N)))
assert(Ext^1(OO_X, N~(>= 0)) != prune truncate(0, Ext^1(truncate(1, S^1), N)))

-- Example 4.2: locally free sheaves and global Ext.
S = ZZ/32003[u,v,w,x,y,z];
I = minors(2, genericSymmetricMatrix(S, u, 3));
X = variety I;
R = ring X;
Omega = cotangentSheaf X;
OmegaDual = dual Omega;
assert(Ext^1(OmegaDual, OO_X^1(>= 0)) == Ext^1(OO_X^1, Omega(>= 0)))

-- Example 4.3: Serre-Grothendieck duality.
S = QQ[v,w,x,y,z];
X = variety ideal(w*x+y*z, w*y+x*z);
R = ring X;
omega = OO_X^{-1};
G = sheaf cokernel genericSymmetricMatrix(R, R_0, 2);
assert(Ext^2(G, omega) == dual HH^0(G))
assert(Ext^1(G, omega) == dual HH^1(G))
assert(Ext^0(G, omega) == dual HH^2(G))
