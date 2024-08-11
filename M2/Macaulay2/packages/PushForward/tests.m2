-- test 0
TEST ///
  kk=ZZ/32003
  R4=kk[a..d]
  R5=kk[a..e]
  R6=kk[a..f]
  M=coker genericMatrix(R6,a,2,3)
  pdim M

  G=map(R6,R5,{a+b+c+d+e+f,b,c,d,e})
  F=map(R5,R4,random(R5^1,R5^{4:-1}))

  P=pushFwd(G,M)
  assert (pdim P==1)

  Q=pushFwd(F,P)
  assert (pdim Q==0)
///

-- test 1
TEST ///
  P3=QQ[a..d]
  M=comodule monomialCurveIdeal(P3,{1,2,3})

  P2=QQ[a,b,c]
  F=map(P3,P2,random(P3^1,P3^{-1,-1,-1}))
  N=pushFwd(F,M)

  assert(hilbertPolynomial M==hilbertPolynomial N)
///

-- test 2
TEST ///
  kk = QQ
  R = kk[x,y]/(x^2-y^3-y^5)
  R' = integralClosure R
  pr = pushFwd map(R',R)
  q = pr_0 / (pr_0)_0
  use R
  assert(ann q==ideal(x,y))
  assert isModuleFinite map(R', R)
///

-- test 3
TEST ///
  kkk=ZZ/23
  kk=frac(kkk[u])
  T=kk[t]
  x=symbol x
  PR=kk[x_0,x_1]
  R=PR/kernel map(T,PR,{t^3-1,t^4-t})
  PS=kk[x_0,x_1,x_2]
  S=PS/kernel map(T,PS,{t^3-1,t^4-t,t^5-t^2})

  rs=map(S,R,{x_0,x_1})
  st=map(T,S,{t^3-1,t^4-t,t^5-t^2})
  assert isModuleFinite rs
  assert isModuleFinite st
  pst=pushFwd st

  MT=pst_0
  k=numgens MT

  un=transpose matrix{{1_S,(k-1):0}}
  MT2=MT**MT

  mtt2=map(MT2,MT,un**id_MT-id_MT**un)
  MMS=kernel mtt2

  r1=trim minimalPresentation kernel pushFwd(rs,mtt2)
  r2=trim minimalPresentation pushFwd(rs,MMS)
  r3=trim (pushFwd rs)_0

  assert(r1==r2)
  assert(flatten entries relations r2 == flatten entries relations r3)
///

-- test 4
TEST ///
  kk=ZZ/3
  T=frac(kk[t])
  A=T[x,y]/(x^2-t*y)

  R=A[p]/(p^3-t^2*x^2)
  S=A[q]/(t^3*(q-1)^6-t^2*x^2)
  f=map(S,R,{t*(q-1)^2})
  assert isModuleFinite f
  pushFwd f

  p=symbol p
  R=A[p_1,p_2]/(p_1^3-t*p_2^2)
  S=A[q]
  f=map(S,R,{t*q^2,t*q^3})
  assert isModuleFinite f
  pushFwd f

  i=ideal(q^2-t*x,q*x*y-t)
  p=pushFwd(f,i/i^3)
  assert(numgens p==2)
///

-- test 5
TEST ///
  kk=QQ
  A=kk[x]
  B=kk[y]/(y^2)
  f=map(B,A,{y})
  assert isModuleFinite f
  pushFwd f
  use B
  d=map(B^1,B^1,matrix{{y^2}})
  assert(pushFwd(f,d)==0)
///

-- test 6
TEST ///
  kk=QQ
  A=kk[t]
  B=kk[x,y]/(x*y)
  use B
  i=ideal(x)
  f=map(B,A,{x})
  assert not isModuleFinite f
  assert(isFreeModule pushFwd(f,module i))
///

-- test 7
TEST ///
  kk=ZZ/101
  n=2

  PA=kk[x_1..x_(2*n)]
  iA=ideal apply(toList(1..n),i->(x_(2*i-1)^i-x_(2*i)^(i+1)))
  A=PA/iA

  PB=kk[y_1..y_(2*n-1)]
  l=apply(toList(1..(2*n-1)),i->(x_i+x_(i+1)))
  g=map(A,PB,l)
  time iB=kernel g;
  B=PB/iB

  f=map(A,B,l)
  assert isModuleFinite f
  assert isModuleFinite g
  time h1=pushFwd g;
  ph1=cokernel promote(relations h1_0,B);
  time h2=pushFwd f;

  assert(ph1==h2_0)
///

-- test 8
TEST ///
  A = QQ
  B = QQ[x]/(x^2)
  N = B^1 ++ (B^1/(x))
  f = map(B,A)
  assert isModuleFinite f
  pN = pushFwd(f,N)
  assert(isFreeModule pN)
  assert(numgens pN == 3)
///

TEST///
  -*
  restart
  *-
  debug needsPackage "PushForward"
  kk = ZZ/101
  A = kk[s]
  B = A[t]
  C = B[u]
  f = map(C,B)
  g = map(C,B,{t})
  assert(isInclusionOfCoefficientRing f)
  assert(isInclusionOfCoefficientRing g)

  kk = ZZ/101
  A = frac (kk[s])
  B = A[t]
  C = B[u]
  f = map(C,B)
  g = map(C,B,{t})
  assert(isInclusionOfCoefficientRing f)
  assert(isInclusionOfCoefficientRing g)
///

TEST///
  -*
  restart
  *-
  debug needsPackage "PushForward"
  s = symbol s; t = symbol t
  kk = ZZ/101
  A = kk[s,t]
  -- note: this ideal is NOT the rational quartic, and in fact has an annihilator over A.
  L = A[symbol b, symbol c, Join => false]/(b*c-s*t, t*b^2-s*c^2, b^3-s*c^2, c^3 - t*b^2)
  isHomogeneous L
  describe L
  basis(L, Variables => L_*)
  inc = map(L, A)
  assert isInclusionOfCoefficientRing inc
  assert isModuleFinite L
  assert isModuleFinite inc
  (M,B,pf) = pushFwd inc
  assert( B*presentation M  == 0)
  assert(numcols B == 5)
///

TEST///
  -*
  restart
  *-
  debug needsPackage "PushForward"
  s = symbol s; t = symbol t
  kk = ZZ/101
  A = kk[s,t]
  L = A[symbol b, symbol c, Join => false]/(b*c-s*t,c^3-b*t^2,s*c^2-b^2*t,b^3-s^2*c)
  isHomogeneous L
  describe L
  basis(L, Variables => L_*)
  inc = map(L, A)
  assert isInclusionOfCoefficientRing inc
  assert isModuleFinite L
  assert isModuleFinite inc
  (M,B,pf) = pushFwd inc -- ok.  this works, but isn't awesome, as it uses a graph ideal.
  assert( B*presentation M  == 0)
  assert(numcols B == 5)
///

TEST///
  -*
  restart
  *-
  debug needsPackage "PushForward"
  s = symbol s; t = symbol t
  kk = ZZ/101
  L = kk[s, symbol b, symbol c, t]/(b*c-s*t, t*b^2-s*c^2, b^3-s*c^2, c^3 - t*b^2)
  A = kk[s,t]
  isHomogeneous L
  inc = map(L, A)
  (M,B,pf) = pushFwd inc
  assert( B * inc presentation M  == 0)
  assert(numcols B == 5)
  pushForward(inc, L^1)
///

TEST///
  -*
  restart
  needsPackage "PushForward"
  *-
  kk = QQ
  A = kk[x]
  R = A[y, Join=> false]/(y^7-x^3-x^2)
  (M,B,pf) = pushFwd map(R,A)
  (M1,B1,pf1) = pushFwd R
  assert(pushFwd(R^3) == pushFwd(map(R,A), R^3))
  assert((M1,B1) == (M,B))
  assert(pushFwd matrix{{y}} == pushFwd(map(R,A),matrix{{y}}))
  assert(isFreeModule M and rank M == 7)
  assert(B == basis(L, Variables => R_*))
  assert( pf(y+x)- matrix {{x}, {1}, {0}, {0}, {0}, {0}, {0}} == 0)
  R' = integralClosure R
  (M,B,pf) = pushFwd map(R',R)
  use R
  assert(M == cokernel(map(R^2,R^{{-6}, {-4}},{{-x^2-x,y^4}, {y^3,-x}})))
  assert(pf w_(2,0) - matrix {{0}, {1}} == 0)
///

TEST ///
  -*
  restart
  needsPackage "PushForward"
  *-
  kk = QQ
  A = kk[x, DegreeRank => 0]
  R = A[y,z, Join => false]
  I = ideal(y^4-x*y-(x^2+1)*z^2, z^4 - (x-1)*y-z^2 - z - y^3)
  B = R/I
  assert isModuleFinite map(B,A)
  (M,g,pf) = pushFwd B
  pushFwd B^1
  pushFwd B^{1}
  fy = pushFwd matrix{{y}}
  fz = pushFwd matrix{{z}}
  assert(fy*fz == pushFwd matrix{{y*z}})
  inc = map(B,A)
  pushFwd(inc, B^1)
  pushFwd(inc, B^{1})
  fy = pushFwd(inc, matrix{{y}})
  fz = pushFwd(inc, matrix{{z}})
  assert(fy*fz == pushFwd(inc, matrix{{y*z}}))

  kk = QQ
  A = kk[x]
  R = A[y,z, Join => false]
  I = ideal(y^4-x*y-(x^2+1)*z^2, z^4 - (x-1)*y-z^2 - z - y^3)
  B = R/I
  assert isModuleFinite map(B,A)
  (M,g,pf) = pushFwd B
  pushFwd B^1
  pushFwd B^{1}
  fy = pushFwd matrix{{y}}
  fz = pushFwd matrix{{z}}
  assert(fy*fz == pushFwd matrix{{y*z}}) -- false
  assert(fy*fz - pushFwd matrix{{y*z}} == 0)
  inc = map(B,A)
  pushFwd(inc, B^1)
  pushFwd(inc, B^{1})
  fy = pushFwd(inc, matrix{{y}})
  fz = pushFwd(inc, matrix{{z}})
  assert(fy*fz == pushFwd(inc, matrix{{y*z}}))

  kk = QQ
  A = kk[x, DegreeRank => 0]
  R = A[y,z]
  I = ideal(y^4-x*y-(x^2+1)*z^2, z^4 - (x-1)*y-z^2 - z - y^3)
  B = R/I
  assert isModuleFinite map(B,A)
  (M,g,pf) = pushFwd B
  pushFwd B^1
  pushFwd B^{1}
  fy = pushFwd matrix{{y}}
  fz = pushFwd matrix{{z}}
  fy*fz == pushFwd matrix{{y*z}}

  kk = QQ
  A = kk[x]
  R = A[y,z]
  I = ideal(y^4-x*y-(x^2+1)*z^2, z^4 - (x-1)*y-z^2 - z - y^3)
  B = R/I
  assert isModuleFinite map(B,A)
  (M,g,pf) = pushFwd B
  pushFwd B^1
  pushFwd B^{{0,1}}
  fy = pushFwd matrix{{y}} -- good
  fz = pushFwd matrix{{z}}
  assert(fy*fz == pushFwd matrix{{y*z}}) -- good
  assert(fy*fz - pushFwd matrix{{y*z}} == 0)
///

TEST ///
  -*
  restart
  needsPackage "PushForward"
  *-
  n = 4
  d = 4
  c = 2
  kk = ZZ/32003;
  S = kk[x_1..x_n];
  I = ideal random(S^1, S^{c:-d});
  R = S/I;
  A = kk[t_1..t_(n-c)];
  phi = map(R, A, random(R^1, R^{n-c:-1}));
  elapsedTime assert isModuleFinite phi
  elapsedTime M1 = pushFwd(phi, R^1)
  elapsedTime M2 = pushForward(phi, R^1);
  assert(M1 == M2)
///

end--

restart
check "PushForward"

-- example bug -----------------------------------
-- DE + MES

///
  restart
  needsPackage "PushForward"

  -- This one works
  kk = ZZ/101
  A = kk[s,t]
  C = A[x,y,z]/(x^2, y^2, z^2)
  phi = map(C,A)
  f = map(C^1, A^4, phi, {{x,s*y,t*y, z}})
  ker f

  -- This one fails, degrees are screwed up.
  kk = ZZ/101
  A = kk[s,t]
  B = frac A
  C = B[x,y,z]/(x^2, y^2, z^2)
  phi = map(C,B)
  f = map(C^1, B^3, phi, {{x,s*y,z}})
  ker f
///

///
  -- Case 1.
  -- ring map is f : A --> B = A[xs]/I, A is a polynomial ring, quotient field, basic field.
///
