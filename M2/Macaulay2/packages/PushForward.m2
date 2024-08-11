---------------------------------------------------------------------------
-- PURPOSE: to compute push forwards of modules and coherent sheaves
--
-- UPDATE HISTORY:
-- - Created Dec 2009 by Claudiu Raicu
-- - Updated Oct 2015 by Karl Schwede
-- - Updated May 2021 by Mike Stillman and David Eisenbud
-- - Updated Aug 2024 by Mahrud Sayrafi
--
-- TODO:
--  finish doc
--  how to interact with pushForward?
--   issues: pushForward seems somewhat faster, in the homogeneous case...
--           also, are these stashed in that case?  (They are not here, yet).
---------------------------------------------------------------------------
newPackage(
    "PushForward",
    Version => "0.8",
    Date => "Aug 11, 2024",
    Authors => {
	{ Name => "Claudiu Raicu",  Email => "craicu@nd.edu",         HomePage => "https://www3.nd.edu/~craicu" },
	{ Name => "David Eisenbud", Email => "de@berkeley.edu",       HomePage => "https://math.berkeley.edu/~de" },
	{ Name => "Mahrud Sayrafi", Email => "mahrud@umn.edu",        HomePage => "https://math.umn.edu/~mahrud/" },
	{ Name => "Mike Stillman",  Email => "mike@math.cornell.edu", HomePage => "https://pi.math.cornell.edu/~mike" }
	},
    DebuggingMode => true,
    AuxiliaryFiles => true,
    Headline => "push forwards of ring maps",
    Keywords => { "Commutative Algebra", "Algebraic Geometry" }
    )

export {
    "isModuleFinite",
    "pushFwd",
    }

-----------------------------------------------------------------------------
-- isModuleFinite
-----------------------------------------------------------------------------

isInclusionOfCoefficientRing = method()
isInclusionOfCoefficientRing RingMap := Boolean => inc -> (
    --checks whether the map is the inclusion of the coefficientRing
    if source inc =!= coefficientRing target inc then return false;
    inc vars source inc == promote (vars source inc, target inc)
    )

isFinite1 = (f) -> (
    A := source f;
    B := target f;
    matB := null;
    mapf := null;
    pols := f.matrix;
    (FA, phiA) := flattenRing A;
    iFA := ideal FA;
    varsA := flatten entries phiA^-1 vars FA;
    RA := try(ring source presentation FA) else FA;
    (FB, phiB) := flattenRing B;
    iFB := ideal FB;
    varsB := flatten entries phiB^-1 vars FB;
    RB := try(ring source presentation FB) else FB;
    m := numgens FA;
    n := numgens FB;
    pols = pols_{0..(m-1)};
    R := try(tensor(RB, RA, Join => false)) else tensor(RB, RA, Join => true);
    xvars := (gens R)_{n..n+m-1};
    yvars := (gens R)_{0..n-1};
    iA := sub(ideal FA,matrix{xvars});
    iB := sub(ideal FB,matrix{yvars});
    iGraph := ideal(matrix{xvars}-sub(pols,matrix{yvars}));
    I := iA+iB+iGraph;
    inI := leadTerm I;
    r := ideal(sub(inI,matrix{yvars | splice{m:0}}));
    for i from 1 to n do
        if ideal(sub(gens r,matrix{{(i-1):0,1_R,(m+n-i):0}}))!=ideal(1_R) then
            return false;
    true
    )

isModuleFinite = method()
isModuleFinite Ring := Boolean => R -> (
    I := ideal leadTerm ideal R;
    ge := flatten select(I_*/support, ell -> #ell == 1);
    set ge === set gens ring I
    )
isModuleFinite RingMap := Boolean => f -> (
    if isInclusionOfCoefficientRing f then
        isModuleFinite target f
    else
        isFinite1 f
    )

-----------------------------------------------------------------------------
-- pushFwd
-----------------------------------------------------------------------------

pushFwd Ring := Sequence => o -> B -> pushFwd(map(B,      coefficientRing B),         o)
pushFwd Module := Module => o -> M -> pushFwd(map(ring M, coefficientRing ring M), M, o)
pushFwd Matrix := Matrix => o -> d -> pushFwd(map(ring d, coefficientRing ring d), d, o)

-- TODO: should this be an internal helper routine?
pushFwd = method(Options => { MinimalGenerators => true })
pushFwd RingMap := Sequence => o -> f -> (
    -- Given ring map f: A --> B
    -- Returns:
    --   pfB is B^1 as an A-module
    --   matB is the set of monomials in B that form a set of generators as an A-module
    --   mapf takes as arg an element of B, and returns a vector in pfB
    (B, A) := (target f, source f);
    (matB, mapfAux) := pushAuxHgs f;
    pushB := makeModule(B^1, f, matB);
    mapf := (b) -> map(pfB, , gens pfB) * mapfAux b;
    (pushB, matB, mapf))

pushFwd(RingMap, Module) := Module => o -> (f, N) -> (
    B := target f;
    aN := ann N;
    C := B/aN;
    bc := map(C, B);
    g := bc * f;
    matB := first pushAuxHgs g;
    M := makeModule(N**C,g,matB);
    if o.MinimalGenerators then prune M else M)

pushFwd(RingMap, Matrix) := Matrix => o -> (f, d) -> (
    (B, A) := (target f, source f);
    src := source d;
    tar := target d;
    --
    C := B / intersect(ann src, ann tar);
    g := map(C, B) * f;
    --
    (matB, mapf) := pushAuxHgs g;
    gR := matB ** cover d;
    --
    psrc := makeModule(src ** C, g, matB);
    ptar := makeModule(tar ** C, g, matB);
    pfd := map(ptar, psrc,
	matrix table(
	    numgens target gR,
	    numgens source gR,
	    (r, c) -> mapf(gR_c_r)));

    if o.MinimalGenerators then prune pfd else pfd)

-- TODO: stash the matB, pf?  Make accessor functions to go to/from gens of R over A, or M to M_A.
-- TODO: given: M = pushFwd N, get the maps from N --> M (i.e. stash it somewhere).
--   also, we want the map going backwards too: given an element of M, lift it to N.

-----------------------------------------------------------------------------
-- makeModule
-----------------------------------------------------------------------------
-- internal function which implements the push forward of a module.
-- input:
--   N     : Module, a module over B
--   f     : RingMap, A --> B
--   matB  : matrix over B, with one row, whose entries form a basis for B over A.
--           in fact, it can be any desired subset of A-generators of B, as well.
-- output:
--   the module N as an A-module.
-- notes:
--   if A is a field, this should be easier?
--   the map mp is basically
--     A^k --> auxN (over B)
--   and its kernel are the A-relations of the elements auxN

makeModule = method()
makeModule(Module, RingMap, Matrix) := (N, f, matB) -> (
     N = trim N;
     auxN := ambient N / image relations N;
     A := source f;
     k := (numgens ambient N) * (numgens source matB);
     --
     mp := if isHomogeneous f
     then try map(auxN, , f, matB ** gens N)
     else map(auxN, A^k, f, matB ** gens N)
     else map(auxN, A^k, f, matB ** gens N);
     --
     ke := kernel mp;
     (super ke) / ke)

-- what if B is an algebra over A (i.e. A is the coefficient ring of B)
-*
    TODO.
    g = gens gb ideal L
    m = lift(matB, ring g)
    coker last coefficients(g, Monomials => m)
*-

-----------------------------------------------------------------------------
-- pushAuxHgs
-----------------------------------------------------------------------------

pushAuxHgs=method()
pushAuxHgs(RingMap):=(f)-> (

    A:=source f;
    B:=target f;

    matB := null;
    mapf := null;

     if isInclusionOfCoefficientRing f then (
     --case when the source of f is the coefficient ring of the target:
	 if not isModuleFinite target f then error "expected a finite map";
	 matB = basis(B, Variables => 0 .. numgens B - 1);
	 mapf = if isHomogeneous f
           then (b) -> (
             (mons,cfs) := coefficients(b,Monomials=>matB);
             lift(cfs, A)
	     )
           else (b) -> (
             (mons,cfs) := coefficients(b,Monomials=>matB);
             -- strip degrees on the target, as otherwise, with differing degrees in A and B,
             -- the degree cannot always be lifted:
             cfs = map(B^(numrows cfs),B^(numcols cfs),cfs);
	     lift(cfs, A)
	     );
         return(matB,mapf)
	 );

     pols:=f.matrix;

     FlatA:=flattenRing A;
     FA:=FlatA_0;
     iFA:=ideal FA;
     varsA:=flatten entries FlatA_1^-1 vars FA;
     RA:=try(ring source presentation FA) else FA;
     FlatB:=flattenRing B;
     FB:=FlatB_0;
     iFB:= ideal FB;
     varsB:=flatten entries FlatB_1^-1 vars FB;
     RB:=try(ring source presentation FB) else FB;
     m:=numgens FA;
     n:=numgens FB;

     pols=pols_{0..(m-1)};

     R := try(tensor(RB, RA, Join => false)) else tensor(RB, RA, Join => true);
     xvars := (gens R)_{n..n+m-1};
     yvars := (gens R)_{0..n-1};
     iA:=sub(ideal FA,matrix{xvars});
     iB:=sub(ideal FB,matrix{yvars});
     iGraph:=ideal(matrix{xvars}-sub(pols,matrix{yvars}));
     I:=iA+iB+iGraph;
     inI:=leadTerm I;

     r:=ideal(sub(inI,matrix{yvars | splice{m:0}}));
     for i from 1 to n do
	if ideal(sub(gens r,matrix{{(i-1):0,1_R,(m+n-i):0}}))!=ideal(1_R) then
     	  error "map is not finite";

     mat:=lift(basis(R/(r+ideal(xvars))),R);
     k:=numgens source mat;
     matB = sub(mat,matrix{varsB|toList(m:0_B)});
     assert(k == numcols matB);
     phi:=map(R,B,matrix{yvars});
     toA:=map(A,R,flatten{n:0_A,varsA});
     mapf = (b)->(
	  (mons,cfs):=coefficients((phi b)%I,Monomials=>mat,Variables=>yvars);
	  toA cfs
	  );
     matB,mapf
     )

-----------------------------------------------------------------------------
-- Tests
-----------------------------------------------------------------------------

load "./PushForward/tests.m2"

-----------------------------------------------------------------------------
-- Documentation
-----------------------------------------------------------------------------

beginDocumentation()

load "./PushForward/docs.m2"

-----------------------------------------------------------------------------
-- Development
-----------------------------------------------------------------------------

end--

restart
uninstallPackage "PushForward"
restart
installPackage "PushForward"
check PushForward
viewHelp PushForward

restart
needsPackage "PushForward"

A = QQ
B = QQ[x]/(x^2)
N = B^1 ++ (B^1/(x))
f = map(B,A)
pushFwd(f,N)
pushFwd f
