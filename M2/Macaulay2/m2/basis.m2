-- Copyright 1995-2002 by Daniel R. Grayson and Michael Stillman
-- Updated 2021 by Mahrud Sayrafi
-* TODO:
 0. hookify, cache
 1. what are (basis, ZZ, List, *) methods for? why only Ring and Ideal?
 2. why isn't (basis, Matrix) implemented?
*-

needs "gb.m2"
needs "max.m2" -- for InfiniteNumber
needs "modules2.m2"
needs "ringmap.m2"

-----------------------------------------------------------------------------
-- Local variables
-----------------------------------------------------------------------------

algorithms := new MutableHashTable from {}

-----------------------------------------------------------------------------
-- Local utilities
-----------------------------------------------------------------------------

inf := t -> if t === infinity then -1 else t

-----------------------------------------------------------------------------
-- helpers for basis
-----------------------------------------------------------------------------

-- the output of this is used in BasisContext
getVarlist = (R, varlist) -> (
    if varlist === null then toList(0 .. numgens R - 1)
    else if instance(varlist, List) then apply(varlist, v ->
        -- TODO: what if R = ZZ?
        if instance(v, R)  then index v else
        if instance(v, ZZ) then v
        else error "expected list of ring variables or variable indices")
    else error "expected list of ring variables or variable indices")

findHeftandVars := (R, varlist, ndegs) -> (
    -- returns (posvars, heftvec)
    -- such that posvars is a subset of varlist
    -- consisting of those vars whose degree is not 0 on the first ndegs slots
    -- and heftvec is an integer vector of length ndegs s.t. heftvec.deg(x) > 0 for each variable x in posvars
    varlist = getVarlist(R, varlist);
    if ndegs == 0 or #varlist == 0 then return (varlist, toList(ndegs : 1));
    if degreeLength R == ndegs and #varlist == numgens R then (
        heftvec := heft R;
        if heftvec =!= null then return (varlist, heftvec));
    --
    zerodeg := toList(ndegs:0);
    posvars := select(varlist, x -> R_x != 0 and take(degree R_x, ndegs) != zerodeg);
    degs := apply(posvars, x -> take(degree R_x, ndegs));
    heftvec = findHeft(degs, DegreeRank => ndegs);
    if heftvec =!= null then (posvars, heftvec)
    else error("heft vector required that is positive on the degrees of the variables " | toString posvars))

-- TODO: can this be done via a tensor or push forward?
-- c.f. https://github.com/Macaulay2/M2/issues/1522
liftBasis = (M, phi, B, offset) -> (
    -- lifts a basis B of M via a ring map phi
    (R, S) := (phi.target, phi.source);
    (n, m) := degreeLength \ (R, S);
    offset  = if offset =!= null then splice offset else toList( n:0 );
    if R === S then return map(M, , B, Degree => offset);
    r := if n === 0 then rank source B else (
        lifter := phi.cache.DegreeLift;
        if not instance(lifter, Function)
        then rank source B else (
            zeroDegree := toList(m:0);
            apply(pack(n, degrees source B),
                deg -> try - lifter(deg - offset) else zeroDegree)));
    map(M, S ^ r, phi, B, Degree => offset))

-----------------------------------------------------------------------------
-- basis
-----------------------------------------------------------------------------

basis = method(TypicalValue => Matrix,
    Options => {
        Strategy   => null,
        SourceRing => null,     -- defaults to ring of the module, but accepts the coefficient ring
        Variables  => null,     -- defaults to the generators of the ring
        Degree     => null,     -- offset the degree of the resulting matrix
        Limit      => infinity, -- upper bound on the number of basis elements to collect
        Truncate   => false     -- TODO: what does this do?
        }
    )

-----------------------------------------------------------------------------

basis Module := opts -> M -> basis(-infinity, infinity, M, opts)
basis Ideal  := opts -> I -> basis(module I, opts)
basis Ring   := opts -> R -> basis(module R, opts)
-- TODO: add? basis Matrix := opts -> m -> basis(-infinity, infinity, m, opts)

-----------------------------------------------------------------------------

basis(List,                           Module) := opts -> (deg,    M) -> basis( deg,   deg,  M, opts)
basis(ZZ,                             Module) := opts -> (deg,    M) -> basis({deg}, {deg}, M, opts)
basis(InfiniteNumber, InfiniteNumber, Module) :=
basis(InfiniteNumber, List,           Module) := opts -> (lo, hi, M) -> basisHelper(opts, lo, hi, M)
basis(InfiniteNumber, ZZ,             Module) := opts -> (lo, hi, M) -> basis( lo,   {hi},  M, opts)
basis(List,           InfiniteNumber, Module) :=
basis(List,           List,           Module) := opts -> (lo, hi, M) -> basisHelper(opts, lo, hi, M)
basis(ZZ,             InfiniteNumber, Module) := opts -> (lo, hi, M) -> basis({lo},   hi,   M, opts)
basis(ZZ,             ZZ,             Module) := opts -> (lo, hi, M) -> basis({lo},  {hi},  M, opts)

-----------------------------------------------------------------------------

basis(List,                           Ideal) :=
basis(ZZ,                             Ideal) := opts -> (deg,    I) -> basis(deg, deg, module I, opts)
basis(InfiniteNumber, InfiniteNumber, Ideal) :=
basis(InfiniteNumber, List,           Ideal) :=
basis(InfiniteNumber, ZZ,             Ideal) :=
basis(List,           InfiniteNumber, Ideal) :=
basis(List,           List,           Ideal) :=
basis(List,           ZZ,             Ideal) :=
basis(ZZ,             InfiniteNumber, Ideal) :=
basis(ZZ,             List,           Ideal) :=
basis(ZZ,             ZZ,             Ideal) := opts -> (lo, hi, I) -> basis(lo,  hi,  module I, opts)

-----------------------------------------------------------------------------

basis(List,                           Ring) :=
basis(ZZ,                             Ring) := opts -> (deg,    R) -> basis(deg, deg, module R, opts)
basis(InfiniteNumber, InfiniteNumber, Ring) :=
basis(InfiniteNumber, List,           Ring) :=
basis(InfiniteNumber, ZZ,             Ring) :=
basis(List,           InfiniteNumber, Ring) :=
basis(List,           List,           Ring) :=
basis(List,           ZZ,             Ring) :=
basis(ZZ,             InfiniteNumber, Ring) :=
basis(ZZ,             List,           Ring) :=
basis(ZZ,             ZZ,             Ring) := opts -> (lo, hi, R) -> basis(lo,  hi,  module R, opts)

-----------------------------------------------------------------------------

basis(List,                           Matrix) :=
basis(ZZ,                             Matrix) := opts -> (deg, M) -> basis(deg, deg, M, opts)
basis(InfiniteNumber, InfiniteNumber, Matrix) :=
basis(InfiniteNumber, List,           Matrix) :=
basis(InfiniteNumber, ZZ,             Matrix) :=
basis(List,           InfiniteNumber, Matrix) :=
basis(List,           List,           Matrix) :=
basis(ZZ,             InfiniteNumber, Matrix) :=
basis(ZZ,             ZZ,             Matrix) := opts -> (lo, hi, M) -> (
    BF := basis(lo, hi, target M, opts);
    BG := basis(lo, hi, source M, opts);
    -- TODO: is this general enough?
    BM := last coefficients(matrix (M * BG), Monomials => BF);
    map(image BF, image BG, BM))

-----------------------------------------------------------------------------

basisHelper = (opts, lo, hi, M) -> (
    R := ring M;
    n := degreeLength R;
    strategy := opts.Strategy;

    -- TODO: check that S is compatible; i.e. there is a map R <- S
    -- perhaps the map should be given as the option instead?
    S := if opts.SourceRing =!= null then opts.SourceRing else R;
    phi := map(R, S);

    if lo === -infinity then lo = {} else
    if lo ===  infinity then error "incongruous lower degree bound: infinity";
    if hi ===  infinity then hi = {} else
    if hi === -infinity then error "incongruous upper degree bound: -infinity";

    -- implement generic degree check
    if #lo != 0  and #lo > n
    or #hi != 0  and #hi > n then error "expected length of degree bound not to exceed that of ring";
    if lo =!= hi and #lo > 1 then error "degree rank > 1 and degree bounds differ";
    if not all(lo, i -> instance(i, ZZ)) then error("expected a list of integers: ", toString lo);
    if not all(hi, i -> instance(i, ZZ)) then error("expected a list of integers: ", toString hi);

    -- e.g., basis(4, 2, QQ[x])
    if #hi == 1 and #lo == 1 and hi - lo < {0}
    then return if S === R then map(M, S^0, {}) else map(M, S^0, phi, {});

    opts = opts ++ {
        Limit      => if opts.Limit == -1 then infinity else opts.Limit
        };

    -- the actual computation of the basis occurs here
    B := runHooks((basis, List, List, Module), (opts, lo, hi, M), Strategy => strategy);

    if B =!= null then liftBasis(M, phi, B, opts.Degree) else if strategy === null
    then error("no applicable strategy for computing bases over ", toString R)
    -- used to be: error "'basis' can't handle this type of ring";
    else error("assumptions for basis strategy ", toString strategy, " are not met"))

-----------------------------------------------------------------------------
-- strategies for basis
-----------------------------------------------------------------------------

basisDefaultStrategy = (opts, lo, hi, M) -> (
    R := ring M;
    A := ultimate(ambient, R); -- is ambient better or coefficientRing?

    -- the assumptions for the default strategy:
    if not (ZZ === A
        or isAffineRing A
        or isPolynomialRing A and isAffineRing coefficientRing A and A.?SkewCommutative
        or isPolynomialRing A and ZZ === coefficientRing A )
    then return null;

    (varlist, heftvec) := findHeftandVars(R, opts.Variables, max(#hi, #lo));

    m := generators gb presentation M;
    log := FunctionApplication { rawBasis, (
            raw m,
            lo, hi,
            heftvec,
            varlist,
            opts.Truncate,
            inf opts.Limit
            )};
    M.cache#"rawBasis log" = Bag {log};
    B := value log;
    B)

-- Note: for now, the strategies must return a RawMatrix
algorithms#(basis, List, List, Module) = new MutableHashTable from {
    Default => basisDefaultStrategy,
    -- TODO: add separate strategies for skew commutative rings, vector spaces, and ZZ-modules
    }

-- Installing hooks for resolution
scan({Default}, strategy ->
    addHook(key := (basis, List, List, Module), algorithms#key#strategy, Strategy => strategy))
