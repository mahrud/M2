-- works automatically: pdim
-- TODO: degree, multidegree
-- TODO: implement the fast powering algorithm
-- TODO: switch the order of the inputs
-- TODO: implement partial caching

importFrom_Core {
    "HilbertContext", "HilbertComputation", "Approximate",
    "isComputationDone", "updateComputation", "adjustComputation"}

getCachedHilbertSeries := M -> if M.cache#?(HilbertContext{}) then M.cache#(HilbertContext{}).Approximate

-- TODO
localPoincareHook = M -> notImplemented "poincare over local rings"
hilbertSamuelPolynomial = opts -> M -> notImplemented "poincare over local rings"

hilbertSamuelSeries = method(Options => options hilbertSeries ++ {cache => null})
hilbertSamuelSeries Module := opts -> M -> (
    -- Output: sum_n H_M(n) * t^n
    if (ord := opts.Order) === infinity then notImplemented "poincare over local rings";
    T := degreesRing ring M;
    coeffs := hilbertSamuelFunction(M, 0, ord-1);
    series := sum(0..ord-1, i -> coeffs_i * T_{i});
    if opts.cache === null then return series else
    container := opts.cache;
    container.Order = ord;
    container.Approximate = series)

-- Computes the Hilbert-Samuel function for modules over local ring, possibly using a parameter ideal.
-- Note:
--   If computing at index n is fast but slows down at n+1, try computing at range (n, n+1).
--   On the other hand, if computing at range (n, n+m) is slow, break up the range.
hilbertSamuelFunction = method()
-- Eisenbud 1995, Chapter 12:
-- Input:  finitely generated (R, m)-module M, integer n
-- Output: H_M(n) := dim_{R/m}( m^n M / m^{n+1} M )
hilbertSamuelFunction(Module, ZZ)     := ZZ   => (M, n)      -> hilbertSamuelFunction(M, {n})
hilbertSamuelFunction(Module, ZZ, ZZ) := List => (M, n0, n1) -> ( --hilbertSamuelFunction(M, {n0}, {n1})
--hilbertSamuelFunction(Module, List)       := ZZ   => (M, L)      -> first hilbertSamuelFunction(M, L, L)
--hilbertSamuelFunction(Module, List, List) := List => (M, L0, L1) -> (
    m := max ring M;
    M = m^n0 * M;
    -*
    wts := heft ring M; -- TODO: must be {1} for this to be correct
    -- series := hilbertSeries(M, Order => 1 + sum(, n1, times));
    series := getCachedHilbertSeries M;
    -- TODO: what's the best algorithm for filling in only the unset ranges?
    if series =!= null then (
	T := monoid ring series;
	(n0', n1') := (n0, n1);
	while coefficient(T_n0, series) != 0 do n0 = n0 + 1;
	while coefficient(T_n1, series) != 0 do n1 = n1 - 1;
	if n1 < n0 then return part(n0', n1', wts, series));
    *-
    for n in n0 .. n1 list (
	if debugLevel >= 1 then printerr("computing HSF_", toString n);
	-- really should be M/mM, but by Nakayama it's the same
	N := localMinimalPresentationHook(M, PruningMap => false);
	if n < n1 then M = m * N;
	numgens N)
    )

hilbertSamuelFunction(Ideal, Module, ZZ)     := ZZ   => (q, M, n) -> first hilbertSamuelFunction(q,M,n,n)
-- Eisenbud 1995, Section 12.1:
-- Input:  parameter ideal q, finitely generated (R,m)-module M, integers n0, n1
-- Output: H_{q, M}(i) := length( q^n M / q^{n+1} M ) for i from n0 to n1
hilbertSamuelFunction(Ideal, Module, ZZ, ZZ) := List => (q, M, n0, n1) -> (
    RP := ring M;
    if class RP =!= LocalRing then error "expected objects over a local ring";
    if ring q =!= RP          then error "expected objects over the same ring";
    if q == max RP            then return hilbertSamuelFunction(M, n0, n1);
    M = localMinimalPresentationHook(M, PruningMap => false);
    M = q^n0 * M;
    for i from n0 to n1 list (
        if debugLevel >= 1 then printerr("computing HSF_", toString i);
        N := localMinimalPresentationHook(M, PruningMap => false);  -- really should be N/mN, but by Nakayama it's the same
        if i < n1 then M = q * N;
        localLengthHook (N/(q * N))
        )
    )

-- TODO: check that it's Artinian first
-- test based on when hilbertSamuelFunction(M, n) == 0?
-- Maybe http://stacks.math.columbia.edu/tag/00IW ?
isFiniteLength = x -> (
    if debugLevel >= 1 then printerr "isFiniteLength is not implemented; assuming the input has finite length"; true)

-- Computes the length of an ideal or module over local rings
-- Note: If computing length is slow, try summing hilbertSamuelFunction for short ranges
-- TODO: implement partial caching
localLengthHook = method(Options => {Strategy => null, Limit => 1000})
localLengthHook Ideal  := ZZ => opts -> I -> localLengthHook(opts, module I)
localLengthHook Module := ZZ => opts -> M -> (
    if not isFiniteLength M then return infinity;
    m := max ring M;
    sum for i to opts.Limit list (
	if M == 0 then break;
        if debugLevel >= 2 then printerr("step ", toString i);
        if opts.Strategy === Hilbert
        then hilbertSamuelFunction(M, i) -- really should be M/mM, but by Nakayama it's the same
        else (
            N := localMinimalPresentationHook(M, PruningMap => false);
            if i < opts.Limit then M = m * N;
            numgens N)
        )
    )

addHook((hilbertFunction, List, Module), Strategy => Default, (L, M) ->
    if instance(ring M, LocalRing) then hilbertSamuelFunction(M, L))

addHook((hilbertSeries, Module), Strategy => Local, (opts, M) ->
    if instance(ring M, LocalRing) then hilbertSamuelSeries(opts, M))

addHook((poincare, Module), Strategy => Local, M ->
    if instance(ring M, LocalRing) then localPoincareHook M)

addHook((hilbertPolynomial, Module), Strategy => Local, (opts, M) ->
    if instance(ring M, LocalRing) then hilbertSamuelPolynomial(opts, M))

end--

restart
errorDepth=2
debug needsPackage "LocalRings"
S = kk[x,y,z]
R = S_(ideal vars S)

elapsedTime hilbertSamuelSeries(module R, Order => 5)
elapsedTime hilbertSeries(module S, Order => 5)
elapsedTime hilbertSeries(module R, Order => 5)


-- TODO
load"~/Projects/M2/Workshops/Workshop-2014-Berkeley/Interpolation/GuessSequence.m2"
see polynomialInterpolation and floaterHormann
