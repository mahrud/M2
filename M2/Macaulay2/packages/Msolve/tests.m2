TEST ///
  -- the msolve hooks must not drop the grading of the target: msolve is handed
  -- bare polynomials, so a basis returned as if the target were S^r shifts every
  -- degree read off `gens gb` -- poincare, hilbertFunction and basis included.
  msolveSetup()
  R = ZZ/3[q]
  N = coker map(R^{-6}, R^{-9}, matrix{{q^3}}) -- generated in degree 6, length 3
  assert(degrees N == {{6}})
  assert(degrees target gens gb presentation N == {{6}})
  assert(degrees source gens gb presentation N == {{9}})
  assert(poincare N == (ring poincare N)_0^6 - (ring poincare N)_0^9)
  assert(apply(10, d -> hilbertFunction(d, N)) == {0,0,0,0,0,0,1,1,1,0})
  assert(degrees source basis N == {{6}, {7}, {8}})
  assert isHomogeneous basis N
///

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
  debug needsPackage "Msolve"
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

TEST ///
  -- msolve's elimination block length is not only an elimination device: it is
  -- a two block degree reverse lexicographic order, and it is Macaulay2's
  -- {GRevLex => a, GRevLex => b} exactly. So a ring with two grevlex blocks is
  -- not msolveApplicable -- msolve cannot be handed the order -- and yet
  -- msolveBlockGB gets the very basis of that order out of it, by computing in
  -- the plain grevlex ring on the same variables.
  debug needsPackage "Msolve"
  kk = ZZ/1073741827
  S = kk[a..e, MonomialOrder => {GRevLex => 2, GRevLex => 3}]
  gensI = {a^2-b*c+d*e, a*b*c-d^2*e+a*e^2, b^3-c*d*e+a^3, c^2*d-a*b*e}

  assert(grevlexBlockLength S === 2)
  assert(grevlexBlockLength(kk[a..e]) === null)
  -- Eliminate n is a weight order in Macaulay2, and a different order
  assert(grevlexBlockLength(kk[a..e, MonomialOrder => Eliminate 2]) === null)
  assert(grevlexBlockLength(kk[a..e, MonomialOrder => {GRevLex => 2, GRevLex => 2}]) === null)

  I = ideal gensI
  assert not msolveApplicable gens I
  ref = sort first entries gens gb I
  G = msolveBlockGB(gens I, msolveDefaultOptions)
  assert instance(G, GroebnerBasis)
  assert(sort first entries gens G == ref)

  -- weights are carried across one block at a time by the same substitution
  W = kk[a..e, Degrees => {2,3,1,4,2},
      MonomialOrder => {GRevLex => {2,3}, GRevLex => {1,4,2}}]
  IW = ideal apply(gensI, f -> sub(f, vars W))
  assert(sort first entries gens msolveBlockGB(gens IW, msolveDefaultOptions)
      == sort first entries gens gb IW)

  -- over a quotient of a block ring the relations go in with the generators and
  -- what in(J) accounts for comes back out, as in msolveGBMatrix
  use S
  R = S/ideal(a^2 - b*c, b^3 - c*d*e)
  M = matrix {{a*b - d*e, c^2*d - a*e^2}}
  assert(sort first entries gens msolveBlockGB(M, msolveDefaultOptions)
      == sort first entries gens gb M)

  -- what msolve cannot do with a block order, and so falls back for: modules,
  -- noncommutative rings, and anything but two blocks
  use S
  assert(msolveBlockGB(matrix {{a,b},{c,d}}, msolveDefaultOptions) === null)
  assert(msolveBlockGB(gens ideal(a^2-b*c, c*d-b^2), msolveDefaultOptions) =!= null)
  -- a plain grevlex ring is msolveForceGB's business, not this one's
  P = kk[a..e]
  assert(msolveBlockGB(substitute(gens I, vars P), msolveDefaultOptions) === null)
  D = ZZ/32003[x,y,dx,dy, WeylAlgebra => {x=>dx, y=>dy},
      MonomialOrder => {GRevLex => 2, GRevLex => 2}]
  assert(msolveBlockGB(matrix {{x*dx - 1, y*dy}}, msolveDefaultOptions) === null)
  E = ZZ/32003[e_1..e_4, SkewCommutative => true,
      MonomialOrder => {GRevLex => 2, GRevLex => 2}]
  assert(msolveBlockGB(matrix {{e_1*e_2, e_3*e_4}}, msolveDefaultOptions) === null)
///

TEST ///
  -- a tower is stored with a block order -- kk[x,y,z][a,b,c] flattens to
  -- kk[a,b,c,x,y,z] with the outer variables their own leading grevlex block,
  -- and that is the order Macaulay2 computes a Groebner basis over the tower in
  -- -- so msolveBlockGB reaches it, the block length being the number of outer
  -- variables and the flattening being carried by flattenRing's two maps.
  debug needsPackage "Msolve"
  kk = ZZ/1073741827
  A = kk[x,y,z]
  B = A[a,b,c]
  assert(grevlexBlockLength ambient first flattenRing B === 3)
  use B

  I = minors_2 matrix {{x,y,z}, {a,b,c}} -- homogeneous in the flattened bigrading
  assert isHomogeneous I
  G = msolveBlockGB(gens I, msolveDefaultOptions)
  assert instance(G, GroebnerBasis)
  assert(ring gens G === B) -- and not the flattening
  assert(sort first entries gens G == sort first entries gens gb I)
  assert isHomogeneous gens G
  assert(degrees target gens G == degrees target gens gb I)
  assert(degrees source gens G == degrees source gens gb I)

  -- a quotient of a tower flattens to a quotient, and its relations go in with
  -- the generators as before
  C = B/ideal(a^2 - b*x)
  M = matrix {{a*b - c*y, b^2*z - a*c*x}}
  assert(sort first entries gens msolveBlockGB(M, msolveDefaultOptions)
      == sort first entries gens gb M)

  -- three levels are three blocks, and msolve has only two
  T = (kk[p])[q][r]
  assert(grevlexBlockLength ambient first flattenRing T === null)
  use T
  assert(msolveBlockGB(matrix {{p*q - r^2, q^3 - p*r}}, msolveDefaultOptions) === null)
///

TEST ///
  -- the gb hook picks a two block grevlex order up, so a quotient ring of such
  -- a ring is built from msolve's basis rather than from Macaulay2's
  debug Core
  debug needsPackage "Msolve"
  kk = ZZ/1073741827
  S = kk[a..e, MonomialOrder => {GRevLex => 2, GRevLex => 3}]
  gensI = {a^2-b*c+d*e, a*b*c-d^2*e+a*e^2, b^3-c*d*e+a^3, c^2*d-a*b*e}
  -- Macaulay2's own answers, taken before the hooks are installed
  ref = sort first entries gens gb ideal gensI
  f = a^3*b + c*d*e^2
  nf = f % gens gb ideal gensI

  msolveSetup {gb}
  I = ideal gensI
  G = gb I
  assert(toString raw G == "declared GB") -- i.e. msolve's, not recomputed
  assert(sort first entries gens G == ref)
  Q = S/I
  assert(lift(sub(f, Q), S) == nf)

  -- other orders still go to Macaulay2's own implementation
  L = kk[x,y,z, MonomialOrder => Lex]
  assert(gens gb ideal(x^2-y*z, y^3-x*z^2, x*y-z^2) != 0)
  T = kk[u_1..u_6, MonomialOrder => {GRevLex => 2, GRevLex => 2, GRevLex => 2}]
  assert(gens gb ideal(u_1^2-u_2*u_3, u_4^2-u_5*u_6, u_1*u_4-u_2*u_5) != 0)
///

TEST ///
  -- the same over a tower, which is where a block order most often comes from:
  -- gb, and so the quotient ring built out of it, goes through msolve
  debug Core
  debug needsPackage "Msolve"
  kk = ZZ/1073741827
  A = kk[x,y,z]
  B = A[a,b,c]
  gensI = first entries gens minors_2 matrix {{x,y,z}, {a,b,c}}
  -- Macaulay2's own answers, taken before the hooks are installed
  ref = sort first entries gens gb ideal gensI
  f = a^2*b*x + c^3*z
  nf = f % gens gb ideal gensI

  msolveSetup {gb}
  I = ideal gensI
  G = gb I
  assert(toString raw G == "declared GB") -- i.e. msolve's, not recomputed
  assert(sort first entries gens G == ref)
  Q = B/I
  assert(lift(sub(f, Q), B) == nf)
  assert(numerator hilbertSeries(Q, Reduce => true)
      == numerator hilbertSeries(B/ideal ref, Reduce => true))
///

TEST ///
  -- the syz hook, and the two options it honors by trimming msolve's answer
  debug needsPackage "Msolve"
  kk = ZZ/1073741827
  R = kk[a..h]
  A = matrix {{a, b, c}, {b, c, d}}
  B = genericMatrix(R, a, 2, 3) | matrix {{a^2}, {b^2}}
  -- inhomogeneous, and so not msolve's: it reads the syzygies off a graded
  -- resolution, which it will not start without a grading
  C = matrix {{a, b+1, c^2}, {b, c, d+a}}
  -- a quotient ring and a block order, neither of which msolve's module F4
  -- takes; note the quotient rebinds a..h, hence the use R
  Q = R/ideal(a^2 - b^2)
  AQ = sub(A, Q)
  P = kk[x,y,z,u,v, MonomialOrder => {GRevLex => 2, GRevLex => 3}]
  AP = matrix {{x, y, z}, {y, z, u}}
  use R

  -- Macaulay2's own answers, taken before the hooks are installed
  refA = syz A
  refB = syz B
  refC = syz C
  refQ = syz AQ
  refP = syz AP
  refnone = syz matrix {{a}, {b}}
  refempty = syz map(R^2, R^0, 0)
  -- rows 1 comes out on the nose; rows 2 only up to a different basis of the
  -- same module, minimal generators being unique only up to that
  refrows0 = syz(B, SyzygyRows => 0)
  refrows1 = syz(B, SyzygyRows => 1)
  refrows2 = syz(B, SyzygyRows => 2)
  -- stopping conditions msolve has no way to impose, so these fall back
  refdeg  = syz(B, DegreeLimit => 3)
  refpair = syz(B, PairLimit => 3)
  refstop = syz(B, StopBeforeComputation => true)

  msolveSetup {syz}
  assert(syz A == refA)
  assert(syz B == refB)
  assert(syz C == refC)
  assert(syz AQ == refQ)
  assert(syz AP == refP)
  assert(syz matrix {{a}, {b}} == refnone)
  assert(syz map(R^2, R^0, 0) == refempty)

  assert(syz(B, SyzygyRows => 1) == refrows1)
  S2 = syz(B, SyzygyRows => 2)
  assert(numrows S2 == 2 and image S2 == image refrows2)
  assert(degrees source S2 == degrees source refrows2)
  assert(syz(B, SyzygyRows => 0) == refrows0)
  -- more rows than there are is not an error, it just asks for all of them
  assert(syz(B, SyzygyRows => 100) == refB)

  -- msolve cannot stop early, so a limit takes that many of the lowest degree
  -- minimal syzygies instead of the ones a stopped computation would hold
  assert(syz(B, SyzygyLimit => 2) == refB_{0,1})
  assert(syz(B, SyzygyLimit => 100) == refB)

  assert(syz(B, DegreeLimit => 3) == refdeg)
  assert(syz(B, PairLimit => 3) == refpair)
  assert(syz(B, StopBeforeComputation => true) == refstop)
///

TEST ///
  -- complex MsolveResolution, and minimizing what comes back.  msolve's
  -- resolution is the nonminimal one, so what has to agree with Macaulay2 is
  -- not its ranks but the minimal Betti numbers underneath them.
  R = ZZ/32003[x,y,z,w]
  I = ideal(x*z-y^2, x*w-y*z, y*w-z^2)
  C = msolveResolution I
  D = complex C
  assert isFree D
  assert(concentration D === (0, length C))
  assert all(1 .. length C, i -> D.dd_i == C.dd_i)
  assert(HH_0 D == coker gens I)
  assert all(1 .. length C - 1, i -> HH_i D == 0)
  assert(betti minimize D == minimalBetti I)

  -- the discriminating one: a rank two free module presented so that msolve
  -- and Macaulay2 disagree about the *nonminimal* resolution, 2,2 against
  -- 2,3,1.  Minimizing has to reconcile them.
  S = ZZ/32003[x,y,z]
  M = map(S^{0,-1}, S^{-2,-2}, {{x^2, y^2}, {z, 0}})
  CM = msolveResolution M
  assert(betti minimize complex CM == minimalBetti coker M)

  -- a frame that runs past the number of variables, where minimizing has
  -- something to cancel at every level
  J = ideal(z, y^2, x^2*y, x^3)
  CJ = msolveResolution J
  assert(length CJ > numgens S)
  assert(betti minimize complex CJ == minimalBetti J)
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
TEST ///
  -- msolve's module F4 has no induced (Schreyer) order, so it must not be
  -- used when the target free module carries one: the result generates the
  -- right submodule but is not a Groebner basis in Macaulay2's order, and
  -- forceGB would declare it one.  Twists of the target are fine, as long as
  -- rawMsolveModuleGB keeps them out of the order it asks msolve for.
  debug needsPackage "Msolve"
  S = ZZ/32003[x,y,z]
  opts = msolveDefaultOptions

  -- constant target degrees: msolve applies, and agrees with Macaulay2
  m = map(S^{-3,-3,-3}, S^{-5,-5,-5}, matrix{{x^2, y*z, 0_S},{z^2, x^2, y^2},{y^2, 0_S, x*z}})
  G = msolveGBMatrix(m, 0, opts)
  assert(G =!= null)
  assert(image leadTerm G == image leadTerm gens gb m)

  -- non-constant target degrees: msolve applies, and the twists are not
  -- allowed to enter the order, so it still agrees with Macaulay2
  mt = map(S^{0,-1}, S^{-2,-3,-2}, matrix{{x^2, y^3, z^2},{x, x*z, y}})
  Gt = msolveGBMatrix(mt, 0, opts)
  assert(Gt =!= null)
  assert(image leadTerm Gt == image leadTerm gens gb mt)

  -- Schreyer order on the target, constant degrees: declined
  F = source schreyerOrder map(S^1, S^{-3,-3,-3}, matrix{{x^3, y^3, z^3}})
  assert(null === msolveGBMatrix(
	  map(F, S^{-5,-5,-5}, matrix{{x^2, y*z, 0_S},{z^2, x^2, y^2},{y^2, 0_S, x*z}}), 0, opts))
///

TEST ///
  -- a Groebner basis over a quotient of a two block grevlex ring: the columns
  -- already accounted for by in(J) are dropped via the quotient ring's own
  -- monomial table, and the result must still match Macaulay2's
  debug needsPackage "Msolve"
  A = ZZ/32003[x,y,z,w, MonomialOrder => {GRevLex => 2, GRevLex => 2}]
  R0 = A/ideal(x*y - z*w)
  gens0 = {x^3 + z^3, y^2*w - x*z, w^4}
  G = gens msolveBlockGB(matrix{gens0}, msolveDefaultOptions)
  GS = gens gb matrix{gens0}
  assert(image G == image GS)
  assert(image leadTerm G == image leadTerm GS)
///

end

restart
needsPackage "Msolve"
check Msolve
