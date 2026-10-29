-- Shared input builders for the msolve benchmark.
-- Every builder is a function of () that constructs its input from scratch, so
-- that no Groebner basis cached on an earlier bmRun's matrix can leak into a
-- later measurement.

bmP = 1073741827  -- prime in msolve's F4SAT window, 2^16 < p < 2^31

-- the ideal of msolve's own input_files/f4sat-31.ms, saturated by z
bmF4satPolys = R -> (
    (w,x,y,z) := (R_0, R_1, R_2, R_3);
    {w^4, x^4, w*y^3 - x*z^3})

bmBaseRing = () -> ZZ/bmP[w,x,y,z]

-- rank 1: the ideal case presented as a module
bmSatRank1 = () -> (
    R := bmBaseRing();
    (image matrix {bmF4satPolys R}, R_3))

-- direct sum of r copies of the rank 1 case: work grows linearly in r while
-- each block stays exactly as hard as the ideal problem
bmSatDirectSum = r -> () -> (
    R := bmBaseRing();
    g := bmF4satPolys R;
    m := matrix table(r, r * #g, (i,j) -> if j // #g == i then g#(j % #g) else 0_R);
    (image m, R_3))

-- rank r with the generators spread across components rather than blocked, so
-- the components genuinely interact
bmSatDense = r -> () -> (
    R := bmBaseRing();
    g := bmF4satPolys R;
    (w,x,y,z) := (R_0, R_1, R_2, R_3);
    cols := flatten apply(#g, k -> apply(r, i -> (
		v := new MutableList from toList(r:0_R);
		v#i = g#k;
		v#((i+1) % r) = (if k == 0 then w else if k == 1 then x else y) * g#((k+1) % #g);
		toList v)));
    (image transpose matrix cols, z))

-- the same shape in higher degree
bmSatDegree = (d, r) -> () -> (
    R := bmBaseRing();
    (w,x,y,z) := (R_0, R_1, R_2, R_3);
    g := {w^d, x^d, w*y^(d-1) - x*z^(d-1)};
    m := matrix table(r, r * #g, (i,j) -> if j // #g == i then g#(j % #g) else 0_R);
    (image m, z))

-- more variables, rank 2
bmSatVars = n -> () -> (
    R := ZZ/bmP[vars(0..n-1)];
    v := gens R;
    g := {v#0^4, v#1^4, v#0 * v#2^3 - v#1 * v#(n-1)^3} |
         apply(toList(3..n-2), i -> v#i^3 * v#0 - v#(i-1)^3 * v#1);
    m := matrix table(2, 2 * #g, (i,j) -> if j // #g == i then g#(j % #g) else 0_R);
    (image m, v#(n-1)))

-- a random submodule, the shape the earlier correctness sweep used
bmSatRandom = (n, r, ngens, d, seed) -> () -> (
    setRandomSeed seed;
    R := ZZ/bmP[vars(0..n-1)];
    m := matrix table(r, ngens, (i,j) -> random(d, R));
    (image m, R_0 + R_1))

-- monomial input: no coefficient growth, so this isolates the combinatorics
bmSatMonomial = (n, r, ngens, d, seed) -> () -> (
    setRandomSeed seed;
    R := ZZ/bmP[vars(0..n-1)];
    mons := flatten entries basis(d, R);
    m := matrix table(r, ngens, (i,j) -> mons#(random(#mons)));
    (image m, R_0))

-- over a quotient ring whose relations are not themselves saturated
bmSatQuotient = r -> () -> (
    A := bmBaseRing();
    (w,x,y,z) := (A_0, A_1, A_2, A_3);
    R := A / ideal(w*z^2, x*y*z);
    g := {sub(w^4, R), sub(x^4, R), sub(w*y^3 - x*z^3, R)};
    m := matrix table(r, r * #g, (i,j) -> if j // #g == i then g#(j % #g) else 0_R);
    (image m, R_3))

-- a subquotient rather than an image
bmSatSubquotient = r -> () -> (
    R := bmBaseRing();
    g := bmF4satPolys R;
    (w,x,y,z) := (R_0, R_1, R_2, R_3);
    m := matrix table(r, r * #g, (i,j) -> if j // #g == i then g#(j % #g) else 0_R);
    rel := matrix table(r, r, (i,j) -> if i == j then w^5 * z else 0_R);
    (subquotient(m, rel), z))

---------------------------------------------------------------------------
-- elimination inputs

bmKatsura = n -> (
    R := ZZ/bmP[x_0..x_n];
    S := i -> if abs i > n then 0_R else x_(abs i);
    L := apply(n, m -> sum(toList(-n..n), i -> S(i) * S(m-i)) - S(m));
    (R, ideal(L | {sum(toList(-n..n), i -> S i) - 1})))

bmCyclic = n -> (
    R := ZZ/bmP[x_0..x_(n-1)];
    L := apply(toList(1..n-1), d ->
	sum(n, i -> product(d, j -> x_((i+j) % n))));
    (R, ideal(L | {product(n, i -> x_i) - 1})))

-- msolve's own input_files/elim-31.ms
bmElimDoc = () -> (
    R := ZZ/bmP[t,w,x,y,z];
    (R, ideal(w^4, x^4, w*y^3 - x*z^3, t*z - 1)))

-- implicitization: the image of a parametrized surface, eliminating parameters
bmImplicitSurface = d -> () -> (
    R := ZZ/bmP[s,t,u,x,y,z];
    (R, ideal(x - s^d, y - t^d, z - s*t*u^(d-2), u^d - s*t)))

bmRandomElim = (n, ngens, d, seed) -> () -> (
    setRandomSeed seed;
    R := ZZ/bmP[vars(0..n-1)];
    (R, ideal apply(ngens, i -> random(d, R))))

---------------------------------------------------------------------------
-- case and method tables, shared by runone.m2 and verify.m2

bmNthread = if getenv "THREADS" == "" then 1 else value getenv "THREADS"

bmSatCases = new HashTable from {
    "rank1"      => bmSatRank1,
    "ds2"        => bmSatDirectSum 2,
    "ds4"        => bmSatDirectSum 4,
    "ds8"        => bmSatDirectSum 8,
    "ds16"       => bmSatDirectSum 16,
    "ds32"       => bmSatDirectSum 32,
    "dense2"     => bmSatDense 2,
    "dense3"     => bmSatDense 3,
    "dense4"     => bmSatDense 4,
    "deg5r2"     => bmSatDegree(5, 2),
    "deg6r2"     => bmSatDegree(6, 2),
    "deg5r4"     => bmSatDegree(5, 4),
    "var5r2"     => bmSatVars 5,
    "var6r2"     => bmSatVars 6,
    "rand4r2"    => bmSatRandom(4, 2, 4, 3, 1001),
    "rand4r3"    => bmSatRandom(4, 3, 5, 3, 1002),
    "rand5r2"    => bmSatRandom(5, 2, 4, 3, 1003),
    "mono4r3"    => bmSatMonomial(4, 3, 8, 5, 2001),
    "mono5r4"    => bmSatMonomial(5, 4, 10, 5, 2002),
    "quot2"      => bmSatQuotient 2,
    "quot4"      => bmSatQuotient 4,
    "subq2"      => bmSatSubquotient 2,
    "subq4"      => bmSatSubquotient 4 }

bmSatMethods = new HashTable from {
    "msolve"         => ((M, f) -> msolveSaturate(M, f, Threads => bmNthread)),
    "msolve-guarded" => ((M, f) -> msolveSaturateVerify(M, f)),
    "m2-iterate"     => ((M, f) -> saturate(M, f, Strategy => Iterate)),
    "m2-default"     => ((M, f) -> saturate(M, f)) }

---------------------------------------------------------------------------
-- each entry is {builder returning (R, I), number of leading vars to eliminate}
bmElimCases = new HashTable from {
    "elimdoc"     => {bmElimDoc, 1},
    "katsura5e2"  => {() -> bmKatsura 5, 2},
    "katsura6e2"  => {() -> bmKatsura 6, 2},
    "katsura7e2"  => {() -> bmKatsura 7, 2},
    "katsura8e2"  => {() -> bmKatsura 8, 2},
    "katsura6e4"  => {() -> bmKatsura 6, 4},
    "katsura7e4"  => {() -> bmKatsura 7, 4},
    "cyclic5e2"   => {() -> bmCyclic 5, 2},
    "cyclic6e2"   => {() -> bmCyclic 6, 2},
    "cyclic7e2"   => {() -> bmCyclic 7, 2},
    "cyclic6e3"   => {() -> bmCyclic 6, 3},
    "surf3"       => {bmImplicitSurface 3, 3},
    "surf4"       => {bmImplicitSurface 4, 3},
    "rand4d3"     => {bmRandomElim(4, 4, 3, 3001), 2},
    "rand5d3"     => {bmRandomElim(5, 5, 3, 3002), 2},
    "rand5d4"     => {bmRandomElim(5, 4, 4, 3003), 2} }

bmElimMethods = new HashTable from {
    -- the path setupMsolve installs for eliminate: it filters the block basis
    -- down to the elimination ideal and returns it in the original ring
    "msolve"     => ((R, I, k) -> msolveEliminationGB(I, take(gens R, k))),
    -- the exported entry point, for reference.  Over ZZ/p with the engine
    -- available this returns the whole block Groebner basis rather than the
    -- elimination ideal, so it is timed but not compared against Macaulay2.
    "msolve-exp" => ((R, I, k) -> msolveEliminate(I, take(gens R, k), Threads => bmNthread)),
    "m2"         => ((R, I, k) -> eliminate(take(gens R, k), I)) }

---------------------------------------------------------------------------
-- module elimination: a free module over a ring where some variables are to be
-- eliminated, given as a matrix whose columns generate the submodule
bmMelimBuild = (n, r, ngens, d, k, seed) -> () -> (
    setRandomSeed seed;
    R := ZZ/bmP[vars(0..n-1)];
    m := matrix table(r, ngens, (i,j) -> random(d, R));
    (R, m, k))

bmMelimCases = new HashTable from {
    "m4r2e1"  => bmMelimBuild(4, 2, 4, 2, 1, 4001),
    "m4r2e2"  => bmMelimBuild(4, 2, 4, 3, 2, 4002),
    "m4r3e2"  => bmMelimBuild(4, 3, 5, 2, 2, 4003),
    "m5r2e2"  => bmMelimBuild(5, 2, 4, 3, 2, 4004),
    "m5r3e2"  => bmMelimBuild(5, 3, 5, 2, 2, 4005),
    "m5r4e3"  => bmMelimBuild(5, 4, 6, 2, 3, 4006),
    "m6r2e3"  => bmMelimBuild(6, 2, 5, 2, 3, 4007),
    "m6r3e3"  => bmMelimBuild(6, 3, 6, 2, 3, 4008) }

-- Macaulay2's own way of doing what msolveElimModuleGB does: a Groebner basis
-- in an order eliminating the first k variables, then the part whose lead term
-- avoids them.  msolveElimModuleGB computes in a two block grevlex ring, which
-- eliminates the same variables, and filters the same way.
bmMelimMethods = new HashTable from {
    "msolve" => ((R, m, k) -> msolveElimModuleGB(m, k, msolveDefaultOptions ++ {Threads => bmNthread})),
    "m2"     => ((R, m, k) -> (
	    n := numgens R;
	    E := (coefficientRing R)(monoid [Variables => n, Degrees => degrees R,
		    MonomialOrder => Eliminate k]);
	    selectInSubring(1, gens gb substitute(m, vars E)))) }

---------------------------------------------------------------------------
-- A canonical form of the answer, computed outside the timed region, so that
-- two processes running different methods on the same case can be compared.

---------------------------------------------------------------------------
-- Parametric homogeneous families
--
-- A case name carries its own parameters, so trying a new size means adding a
-- line to a spec file rather than an entry here: hsat-n5-r3-k6-d3-s7 is five
-- variables, three components, six generators of degree three, seed seven.
--
-- A matrix whose entries are all forms of one degree is homogeneous with the
-- target free module in degree zero, which is the shape every family here
-- uses, so every module below is graded and so is the element saturated by.

bmParams = s -> new HashTable from apply(drop(separate("-", s), 1),
    p -> substring(p, 0, 1) => value substring(p, 1))

bmHomogMatrix = (R, r, k, d, mono) -> (
    if not mono then matrix table(r, k, (i, j) -> random(d, R))
    else (
	B := flatten entries basis(d, R);
	matrix table(r, k, (i, j) -> B#(random(#B)))))

bmHomogRing = p -> ZZ/bmP[vars(0 .. p#"n" - 1)]

-- random forms, saturated by a linear form
bmHsat = nm -> () -> (
    p := bmParams nm;
    setRandomSeed p#"s";
    R := bmHomogRing p;
    (image bmHomogMatrix(R, p#"r", p#"k", p#"d", false), R_0 + R_1))

-- random monomials, saturated by a variable: no coefficient growth, so this
-- isolates the combinatorics, and it was the hardest shape in the first sweep
bmHmon = nm -> () -> (
    p := bmParams nm;
    setRandomSeed p#"s";
    R := bmHomogRing p;
    (image bmHomogMatrix(R, p#"r", p#"k", p#"d", true), R_0))

-- homogeneous module elimination: e leading variables to eliminate
bmHmel = nm -> () -> (
    p := bmParams nm;
    setRandomSeed p#"s";
    R := bmHomogRing p;
    (R, bmHomogMatrix(R, p#"r", p#"k", p#"d", false), p#"e"))

bmLookupSat = nm -> (
    if bmSatCases#?nm then bmSatCases#nm
    else if match("^hsat-", nm) then bmHsat nm
    else if match("^hmon-", nm) then bmHmon nm
    else error("unknown saturation case " | nm))

bmLookupMelim = nm -> (
    if bmMelimCases#?nm then bmMelimCases#nm
    else if match("^hmel-", nm) then bmHmel nm
    else error("unknown module elimination case " | nm))
