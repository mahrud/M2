TEST ///
  -- c.f. https://github.com/algebraic-solving/msolve/issues/165
  K = ZZ/65537 -- > 2^16
  R = K[x_(0,0), x_(0,1)]
  I = ideal(2*x_(0,0)+3*x_(0,1))
  assert(I == ideal msolveSaturate(I, x_(0,0)))
  assert(I == saturate(I, x_(0,0), Strategy => Msolve))

  K = ZZ/1073741827 -- > 2^30
  R = K[x_(0,0), x_(0,1)]
  I = ideal(2*x_(0,0)+3*x_(0,1))
  assert(I == ideal msolveSaturate(I, x_(0,0)))
  assert(I == saturate(I, x_(0,0), Strategy => Msolve))

  -- F4SAT doesn't work with small primes
  K = ZZ/32771 -- < 2^16
  R = K[x_(0,0), x_(0,1)]
  I = ideal(2*x_(0,0)+3*x_(0,1))
  assert(try (ideal msolveSaturate(I, x_(0,0));         false) else true)
  assert(try (saturate(I, x_(0,0), Strategy => Msolve); false) else true)
  assert(I == saturate(I, x_(0,0)))
///

TEST ///
  A = ZZ/1073741827[x,y,z]
  B = A[u,v,w]
  I = minors_2 matrix {{x,y,z}, {u,v,w}}
  -- TODO: loses homogeneity
  assert(I == ideal msolveGB I)
  assert(I == ideal msolveSaturate(I, B_0))
  assert(I == ideal msolveSaturate(I, B_3))
  -- TODO: eliminate doesn't work over tower rings
  msolveEliminate(I, B_0) -- == eliminate(I, B_0)
  -- TODO: msolveEliminate can't eliminate variables in the base ring:
  -- msolveEliminate(I, B_3)
  -- TODO: msolveLeadMonomials doesn't work for tower rings yet
  -- msolveLeadMonomials I
  --
  I = minors_2 matrix {{x,y^2,z^3}, {u^4,v^5,w^6}}
  assert(I == ideal msolveGB I)
///

TEST ///
  R = QQ[x,a,b,c,d]
  f = 7*x^2+a*x+b
  g = 2*x^2+c*x+d
  I = eliminate(x, ideal(f,g))
  J = msolveEliminate(x, ideal(f,g))
  assert(sub(I, ring J) == J)

  R = ZZ/11[x,a,b,c,d]
  f = 7*x^2+a*x+b
  g = 2*x^2+c*x+d
  I = eliminate(x, ideal(f,g))
  J = msolveEliminate(x, ideal(f,g))
  assert(sub(I, ring J) == ideal selectInSubring(1, gens J))
///

TEST ///
  R = ZZ/1073741827[z_1..z_3]
  I = ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
  assert(leadTerm msolveGB I == matrix{{z_1*z_2, z_1^2, z_1*z_3^2, z_2^2*z_3^2}})

  R = QQ[z_1..z_3]
  I = ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
  assert(ideal msolveGB I == ideal groebnerBasis I)
  assert(leadTerm msolveGB I == matrix{{7*z_1*z_2, 8*z_1^2, 56*z_1*z_3^2, 235*z_2^2*z_3^2}})
///

TEST ///
  -- msolveGBMatrix over a quotient ring R = S/J, fed a matrix with more than
  -- one row: rels = presentation R must be tensored with id_(target m0) so
  -- its row count matches m0's, and the in(J) filter afterward must inspect
  -- the whole lead term column, not just its first row, or generators whose
  -- lead term is nonzero only past row 1 are silently dropped.
  S = ZZ/32003[x,y,z, MonomialOrder => GRevLex]
  J = ideal(x^2 - y*z, y^3 - x*z^2)
  R = S/J
  M = matrix {{x, y}, {y, z}, {z, x^2}}
  G = msolveGB M
  assert(image G == image M)
  assert(image G == image gens gb M)
///

TEST ///
  -- rawMsolveGB dispatches a multi-row matrix to msolve's module F4, using
  -- Macaulay2's default module order: term over position up.
  debug Core
  needsPackage "Msolve"
  S = ZZ/32003[x,y,z,w]
  A = matrix {{x,y,z}, {y,z,w}}
  I = minors_2 A
  assert(msolveGB I == gens gb I)
  assert(msolveGB A == gens gb A)
  assert(msolveSyzygy A === syz A)
  assert(unpackMsolveBetti rawMsolveMinimalBetti(raw A, 5, 1, 0) === minimalBetti coker A)
  assert(rawMsolvePoincare(raw A, 3, 0) === raw poincare coker A)

  -- technically the output is a gb but not mingens
  -- does M2 do extra work to get the mingens? might be worth it
  S = ZZ/32003[x,y,z,w]
  A = matrix {{x,y,z}, {y,z,w}}
  I = minors_2 A
  msolveSyzygy gens I
  syz gb(gens I, Syzygies => true)
  res(coker gens ideal I_*, Strategy => Nonminimal)

  gens gb syz gens truncate(1, S)
  msolveSyzygy gens truncate(1, S) -- FIXME

  R = quotient I
  syz gens truncate(2, R)
  msolveSyzygy gens truncate(2, R) -- FIXME

  -- Claude profile this example
  restart
  debug Core
  debug needsPackage "Msolve"
  needsPackage "NormalToricVarieties"
  X = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/101)
  S = ring X
  M = truncate({4,4,4}, S^2);
  f = raw presentation M;
  elapsedTime assert(unpackMsolveBetti rawMsolveMinimalBetti(f, 3, 1, 0) == betti res M) -- 1.5s
  elapsedTime assert(rawMsolvePoincare(f, 1, 0) === raw poincare M) -- 1.4s

  -- Claude profile this advanced example
  m = 4 -- >= 4
  d = 2*m-1
  N = ZZ^d
  A = map(N, ZZ^1, 0) | N_{0..d-2} | (sum(d-1, i -> (m-1) * N_{i}) + m*N_{d-1})
  P = 3 * convexHull A
  D = toricDivisor(P, CoefficientRing => ZZ/101)
  X = variety D
  S = ring X
  isVeryAmple D   -- isVeryAmple P
  isNormal P
  isSimplicial X  -- isSimplicial P
  classGroup X
  effGenerators S -- matrix transpose degrees S
  I = ideal monomials D
  -- elapsedTime C = res I
  -- elapsedTime poincare I;
  -- elapsedTime minimalBetti I; -- ??
  errorDepth = 2
  elapsedTime unpackMsolveBetti rawMsolveMinimalBetti(raw gens I, 8, 1, 2);
  elapsedTime assert(rawMsolvePoincare(raw gens I, 1, 2) === raw poincare I)
///
	      
TEST ///
  R = QQ[x,y];
  I = ideal ((x-3)*(x^2+1),y-1);
  assert(msolveRealSolutions I == {{3.0, 1.0}})
  assert(msolveRealSolutions I === msolveRealSolutions(I, QQi))
  scan({QQ, QQi, RR, RR_53, RRi, RR_53},
      F -> assert({{3.0, 1.0}} == msolveRealSolutions(I, F)))
  assert(precision first first msolveRealSolutions I           == infinity)
  assert(precision first first msolveRealSolutions(I, RR)      == defaultPrecision)
  assert(precision first first msolveRealSolutions(I, RR_20)   == defaultPrecision)
  assert(precision first first msolveRealSolutions(I, RR_100)  == 100)
  assert(precision first first msolveRealSolutions(I, RRi)     == defaultPrecision)
  assert(precision first first msolveRealSolutions(I, RRi_20)  == defaultPrecision)
  assert(precision first first msolveRealSolutions(I, RRi_100) == 100)
///

TEST ///
  S = ZZ/1073741827[t12,t13,t14,t23,t24,t34];
  I = ideal(
      (t13*t12-t23)*(1-t14)+(t14*t12-t24)*(1-t13) - (t12+1)*(1-t13)*(1-t14),
      (t23*t12-t13)*(1-t24)+(t24*t12-t14)*(1-t23) - (t12+1)*(1-t23)*(1-t24),
      (t14*t34-t13)*(1-t24)+(t24*t34-t23)*(1-t14) - (t34+1)*(1-t14)*(1-t24));
  sat = (1-t24)
  elapsedTime J1 = saturate(ideal I_*, sat, Strategy => Eliminate);
  elapsedTime J2 = ideal msolveSaturate(ideal I_*, sat);
  assert(J1 == J2)
///

TEST ///
  R = QQ[x,a,b,c,d];
  f = x^2+a*x+b;
  g = x^2+c*x+d;
  I = ideal(f, g)
  eM2 = eliminate(x, I);
  eMsolve = msolveEliminate(x, I);
  assert(eM2 == sub(eMsolve, ring eM2))
///

TEST ///
R=QQ[x,y,z];
I = ideal {(x^3-z)*(x^2-y)};
S = ideal {z-2,y-1};
rur = msolveRUR(I+S);
p = matrix{{1}};
assert(apply(rur#"numerator", n->sub(-n,p)) / sub(rur#"denominator",p)=={1,1,2});
p = matrix{{-1}};
assert(apply(rur#"numerator", n->sub(-n,p)) / sub(rur#"denominator",p)=={-1,1,2});
///

TEST ///
  debugLevel=1
  R = QQ[x..z,t]
  K = ideal(x^6+y^6+x^4*z*t+z^3,36*x^5+60*y^5+24*x^3*z*t,
      -84*x^5+10*x^4*t-56*x^3*z*t+30*z^2,-84*y^5-6*x^4*t-18*z^2,
      48*x^5+10*x^4*z+32*x^3*z*t,48*y^5-6*x^4*z,14*x^4*z+8*x^4*t+24*z^2)
  errorDepth=2
  W1 = msolveEliminate(R_0, K, Verbosity => 1)
  W2 = eliminate(R_0, K)
  assert(sub(W1, R) == W2)
///

TEST ///
  -- genus 6, 14th semigroup
  msolveSetup()
  kk = ZZ/32003
  R = kk[x_0..x_7]/(x_5^2-x_4*x_6,x_4*x_5-x_3*x_7,x_3*x_5-x_2*x_6,x_4^2-x_2*x_7,x_3*x_4-x_2*x_5,
      x_2*x_4-x_1*x_7,x_3^2-x_1*x_6,x_2*x_3-x_1*x_5,x_1*x_3-x_0*x_6,x_0*x_3-5*x_6^2+5*x_5*x_7,
      x_2^2-x_1*x_4,x_1*x_2-x_0*x_5,x_0*x_2-5*x_5*x_6+5*x_4*x_7,x_1^2-5*x_6^2+5*x_5*x_7,
      x_0*x_1-5*x_3*x_6+5*x_2*x_7,x_0^2*x_6-5*x_1*x_6^2+5*x_1*x_5*x_7,
      x_0^2*x_5-5*x_1*x_5*x_6+5*x_1*x_4*x_7,x_0^2*x_4-5*x_1*x_4*x_6+5*x_0*x_7^2)
  I = ideal(-x_0,-x_1,-x_2,-x_3,-x_4,-x_5,-x_6)
  L = map(R^1, R^{5:{-10}, 2:{-20}}, {{
              x_6^10+4*x_5*x_6^8*x_7-x_4*x_6^7*x_7^2-2*x_3*x_6^5*x_7^4+2*x_2*x_6^4*x_7^5+x_2*x_5*x_6^2*x_7^6-2*x_0^2*x_7^8-5*x_1*x_6*x_7^8,
              x_5*x_6^9+4*x_4*x_6^8*x_7+x_0^7*x_7^3-2*x_3*x_6^6*x_7^3+3*x_2*x_6^5*x_7^4+3*x_2*x_5*x_6^3*x_7^5+2*x_1*x_5*x_7^8,
              x_3*x_6^9+4*x_2*x_6^8*x_7-x_2*x_5*x_6^6*x_7^2-5*x_0^6*x_7^4-x_1*x_6^5*x_7^4-2*x_1*x_5*x_6^3*x_7^5-4*x_1*x_4*x_6^2*x_7^6+3*x_0*x_6*x_7^8,
              x_0*x_4*x_6^8-x_0^8*x_7^2-3*x_6^8*x_7^2-2*x_5*x_6^6*x_7^3-x_4*x_6^5*x_7^4-2*x_3*x_6^3*x_7^6+3*x_2*x_6^2*x_7^7+5*x_2*x_5*x_7^8,
              x_0^10, x_0^19*x_7, x_0*x_5*x_6^18-2*x_0*x_4*x_6^17*x_7+5*x_6^17*x_7^3-5*x_5*x_6^15*x_7^4-4*x_0^5*x_7^15}})
  a = gcd flatten degrees source L
  degs = degrees source L // a
  T = kk(monoid[ #degs, VariableBaseName => getSymbol "s", Degrees => degs ]);
  f = map(R, T, L, DegreeMap => d -> d * a)
  assert isHomogeneous f
  -- allowableThreads = 8
  gbTrace = 2
  -- you should see msolve's internal log
  elapsedTime J = ker f; -- 23s
  assert(numgens J == 972);
///

/// 
kk = QQ
R1 = kk[a..f, MonomialSize=>8];
setRandomSeed 42
J1 = ideal random(R1^1, R1^{-2,-2,-3,-3}, Height=>100);
elapsedTime gbC = flatten entries gens (G = gb(ideal J1_*));
gbMsolve = flatten entries msolveGB ideal J1_*;
assert(gbC == gbMsolve)
///
end

restart
needsPackage "Msolve"
check Msolve
