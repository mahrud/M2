-- TODO: add assertions
TEST /// -- test solveFrobeniusIdeal
  S = QQ[x_1..x_5]
  W = makeWeylAlgebra S
  T = first createThetaRing W
  -- see SST Example 2.3.16
  w = {1,1,1,1,1}
  J = ideal(T_0+T_1+T_2+T_3+T_4, T_0+T_1-T_3, T_1+T_2-T_3, T_0*T_2, T_1*T_3)
  F = solveFrobeniusIdeal J
  g = map(W, T, apply(5, i -> W_i*W_(i+5)))
  cssLeadTerms(g J, w)
  -- FIXME
  --truncatedCanonicalSeries(g J, w, 5)
///

TEST ///
  S = QQ[x]
  w = {1}
  W = makeWeylAlgebra S;

  I = ideal(x*dx*(x*dx-3)-x*(x*dx+101)*(x*dx+13))
  assert(nilssonSupport(I,w,3) == {{0},{1},{2},{3}})
  assert(toString cssLeadTerms(I, w) == "{1, X_0^3}")
  (G, sols) = truncatedCanonicalSeries(I, w, 4)
  residues = table(G, sols, applyNilssonOperator)
  assert all(flatten residues, f -> all(exponents f, e -> nilssonWeight_w e > 4))

  I = ideal(x*dx*(x*dx-3) - x*(x*dx+10)*(x*dx+20))
  -- this used to crash when k = 3; see https://github.com/Macaulay2/M2/issues/2831
  (G, sols) = truncatedCanonicalSeries(I, w, 3)
  residues = table(G, sols, applyNilssonOperator)
  assert all(flatten residues, f -> all(exponents f, e -> nilssonWeight_w e > 3))
  (G, sols) = truncatedCanonicalSeries(I, w, 4)
  residues = table(G, sols, applyNilssonOperator)
  assert all(flatten residues, f -> all(exponents f, e -> nilssonWeight_w e > 4))
///

TEST ///
  debug needsPackage "HolonomicSystems" -- for nonpositiveWeightGens
  -- simple version first
  W = makeWeylAlgebra(QQ[x_0,x_1]); w = {-10,-17}
  netList(I = ideal(x_0*dx_0^2 - x_1*dx_1^2 + dx_0 - dx_1, x_0*dx_0 + x_1*dx_1 + 1))_*
  netList(J = ideal nonpositiveWeightGens(I, w))_*
  assert(toString cssLeadTerms(I, w) == "{Xinv_1, -Xinv_1*logX_0+Xinv_1*logX_1}")
  elapsedTime (G, sols) = truncatedCanonicalSeries(I, w, 10); -- <1s
  residues = table(G, sols, applyNilssonOperator)
  assert all(flatten residues, f -> all(exponents f, e -> nilssonWeight_w e > 30))

  -- Lizzie Pratt's example
  W = makeWeylAlgebra(QQ[x_1,x_2,x_3]); w = {-1,0,1}
  netList(I = ideal(
      x_1*dx_1^2 - x_3*dx_3^2 + dx_1 - dx_3,
      x_2*dx_2^2 - x_3*dx_3^2 + dx_2 - dx_3,
      x_1*dx_1 + x_2*dx_2 + x_3*dx_3 + 1))_*
  netList(J = ideal nonpositiveWeightGens(I, w))_*
  assert(sort nilssonSupport(J, w, 2) == {{-2,2,0},{-1,0,1},{-1,1,0},{0,0,0}})
  assert(toString cssLeadTerms(I, w) == "{Xinv_0, -Xinv_0*logX_0+Xinv_0*logX_2, -Xinv_0*logX_0+Xinv_0*logX_1, Xinv_0*logX_0^2-Xinv_0*logX_1*logX_0-Xinv_0*logX_2*logX_0+Xinv_0*logX_2*logX_1}")
  elapsedTime (G, sols) = truncatedCanonicalSeries(I, w, 2) -- ~20s
  residues = table(G, sols, applyNilssonOperator)
  assert all(flatten residues, f -> all(exponents f, e -> nilssonWeight_w e > 2))
///

TEST ///
  W = makeWA(QQ[x_1,x_2])
  assert all({0,1}, b -> isTorusFixed ideal(
          x_1*dx_1*(x_1*dx_1+b), x_1*dx_1*(x_2*dx_2+b),
          x_2*dx_2*(x_1*dx_1+b), x_2*dx_2*(x_2*dx_2+b)))
  assert all({0,1}, b -> isTorusFixed ideal(
          x_1*dx_1*x_2*dx_2-1, x_1*dx_1-b*x_2*dx_2-b+1))
///

TEST ///
  W = makeWeylAlgebra(QQ[x,y])
  vars W
  thetax = x*dx
  thetay = y*dy
  P1 = thetax^2*(thetax-2)-x*(thetax+thetay+1)*(thetax+2)*(thetax+3)
  P2 = thetay^2*(thetay-3)-y*(thetax+thetay+1)*(thetay+2)*(thetay+3)
  I = ideal(P1,P2)
  w = {1,11}
  inw(I,flatten{-w|w})
  -- elapsedTime (G, sols) = truncatedCanonicalSeries(I, w, 2)
  -- residues = table(G, sols, applyNilssonOperator)
  -- assert all(flatten residues, f -> all(exponents f, e -> nilssonWeight_w e > 2))

  S = QQ[t_1,t_2]
  distraction(I,S)
  cssExptsMult(I,w)
  --{{4, {0, 0}}, {2, {2, 0}}, {2, {0, 3}}, {1, {2, 3}}}
  --matches SST Ex 2.5.13
///

TEST ///
  A = matrix{{1,1,1,1,1},{1,1,0,-1,0},{0,1,1,-1,0}}
  beta = {1,0,0}
  I = gkz(A,beta)
  w = {1,1,1,1,0}
  S = QQ[t_1..t_5]
  isTorusFixed I --false
  J = inw(I,flatten{-w|w}) 
  isTorusFixed J --true
  distraction(J,S) == ideal(t_1 +t_2 +t_3+t_4 +t_5 -1, t_1 +t_2 -t_4, t_2 +t_3 -t_4, t_1*t_3, t_2*t_4)
  cssExptsMult(I,w) --{{4, {0, 0, 0, 0, 1}}}
  --matches Ex 2.6.4
///

-- FIXME: nilssonSupport fails here, because truncating the Nilsson cone at
-- weight k leaves an unbounded polyhedron and latticePoints gives up with
-- "Something went wrong, vertex with negative height."
-- TODO: do SST eq. (1.22); see SST pp. 26
-*
  A = matrix{{1,0,0,-1},{0,1,0,1},{0,0,1,1}}
  beta = {1,0,0}
  I = gkz(A,beta)
  w = {1,1,1,1} -- NOTE: this used to be {1,1,1,1,0}, which is not a weight for I
  nilssonSupport(I,w)
  nilssonSupport(I,w,3)
  cssLeadTerms(I, w)
  (G, sols) = truncatedCanonicalSeries(I, w, 4);
*-

end--
restart
needsPackage "HolonomicSystems"
check HolonomicSystems
viewHelp HolonomicSystems

