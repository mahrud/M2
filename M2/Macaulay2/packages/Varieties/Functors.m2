-- Conventions:
-- - R is the coordinate ring of X = Proj R or Spec R
-- - S is a polynomial ring (= symmetric algebra)
--   i.e. coordinate ring of the ambient space PP^n or AA^n.
-- -

-----------------------------------------------------------------------------
-- helpers for sheaf cohomology
-----------------------------------------------------------------------------

-- FIXME: a kludge to allow the package to load with OldChainComplexes
if HomologicalAlgebraPackage =!= "Complexes" then (
    protect Concentration;
    freeResolution = res;
    koszulComplex = x -> notImplemented();
    liftComplex = x -> notImplemented();
    )

-- TODO: should this also check that the variety is finite type over the field?
checkVariety := (X, F) -> (
    if not X === variety F     then error "expected coherent sheaf over the same variety";
    if not isAffineRing ring X then error "expected a variety defined over a field";
    )

degreeList := M -> (
    -- gives the exponents of the numerator of reduced Hilbert series of M
    if dim M > 0 then error "expected module of finite length";
    H := poincare M; -- HS(M) = poincare M / product_i (1 - T^(d_i)),
    T := ring H;
    -- this is faster than calling reduceHilbert
    w := product(degrees ring M, deg -> 1 - T_deg);
    first \ exponents(H // w))

-- quotienting by the local cohomology H_m^0(M) to "saturate" M
-- TODO: use irrelevant ideal here
killLocalH0 = M -> M.cache.TorsionFree ??= if (H0 := saturate(0*M)) == 0 then M else M / H0
-- We mainly remember M.cache.TorsionFree (as a quotient module of M) for use in SheafMaps.m2.

protect SaturatedLift
protect TorsionFreeLift
protect TorsionFreeMap

-- Given a sheaf F associated to a graded module M over a positively graded algebra R,
-- return the map from M to a possibly simpler R-module N that represents the same sheaf F.
-- We always simplify at least to M/M_tors, and if an even "better" module has been cached, we return that.
-- The module N can be obtained by "target torsionFreeMap F".
torsionFreeMap = F -> (
    M := module F;
    R := ring M;
    if F.cache.?SaturationMap  then F.cache.SaturationMap  else
    if F.cache.?TorsionFreeMap then F.cache.TorsionFreeMap else (
	N0 := killLocalH0 M;
	N := minimalPresentation N0;
	F.cache.TorsionFreeLift = minimalPresentation liftModule N;
	F.cache.TorsionFreeMap = inverse N.cache.pruningMap * inducedMap(N0, M)
	)
    )

torsionFreeLift = F -> (
    M := module F;
    if F.cache.?SaturatedLift   then F.cache.SaturatedLift   else
    if F.cache.?TorsionFreeLift then F.cache.TorsionFreeLift else (
	N0 := killLocalH0 M;
	N := minimalPresentation N0;
	F.cache.TorsionFreeMap = inverse N.cache.pruningMap * inducedMap(N0, M);
	F.cache.TorsionFreeLift = minimalPresentation liftModule N)
    )

-- This function "twistedGlobalSectionsModule" returns a sequence (b0, b1, M0) as output,
-- meaning that M0 is a module that maps to HH^0(X, F(*)),
-- the map is an isomorphism in degrees at least b0, and it is surjective in degrees at least b1.
-- (These bounds need not be optimal.) Here F can be a coherent sheaf on a closed subspace
-- of a weighted projective space.
--
-- TODO: add tests:
-- - global sections of sheafHom are Hom
-- TODO: implement for multigraded ring using emsbound
-- TODO: this can change F.module to the result!
-- NOTE: this may have elements in degrees below bound as well, is that a bug?
twistedGlobalSectionsModule = (F, bound) -> (
    -- compute global sections module Gamma_(d >= bound)(X, F(d)), for a CoherentSheaf F.
    -- But if HH^0(X, F(*)) is bounded below, then return the complete answer, regardless of "bound".
    -- FIXME: this line, as opposed to
    --  cokernel presentation module F
    -- breaks the test added in 975d780470.
    -- However, we need to keep the information
    -- cached in M, for instance if M is a Hom module.
    M := module F;
    R := assertWeightedZZGraded ring M;
    quot := torsionFreeMap F; -- The map from M to a simplified R-module N (at least simplified to M/M_tors).
    N := target quot;
    N' := torsionFreeLift F; -- This is N as an S-module.
    S := ring N'; -- X is a closed subvariety in Proj S
    n := numgens S; -- 1 + dim Proj S
    w := S.cache.Dualizing ??= S^{- sum degrees S};
    -- We fix the dualizing module w, as a graded S-module.
    -- As a result, Ext^j_S(N', w), in case another function has computed that module earlier.
    --
    -- Note: bound=infinity may mean that HH^1_m(N) = 0, i.e., N is saturated,
    -- or just that the user does not want to search for any global sections besides those in N.
    -- TODO: what would pdim N' < n-1, hence E1 = 0, imply?
    complete := false;
    p := 0;
    if bound < infinity then (
	p = if pdim N' < n-1 then (
	    complete = true;
	    0)
	else (
	    E1 := Ext^(n-1)(N', w); -- This is the graded dual of the local cohomology HH^1_m(N').
	    -- For example, if N' = M/M_tors, then this is isomorphic to HH^1_m(M).
	    if dim E1 <= 0 then ( -- 0-module or 0-dim module (i.e., finite length))
		complete = true;
		degList := degreeList E1;
		1 + max degList - min degList
		)
	    else 1 - first min exponents poincare E1 - bound
	    -- Note that "first min degrees E1" would give the smallest degree of a generator for E1,
	    -- even if that generator is killed by the relations. Instead, we use poincare E1 to find the smallest degree
	    -- in which E1 is not zero. For example, if poincare E1 = 5T^(-1)+T^2, then exponents poincare E_1 = {{-1},{2}},
	    -- hence the need for the "first".
	    )
	);
    -- this can only happen if bound=-infinity, that is, from calling HH^0(F(*)) = HH^0(F(>=(-infinity))
    if p === infinity then error "the global sections module is not finitely generated";
    -- caching these to be used later in prune SheafMap
    F.cache.TorsionFree = M.cache.TorsionFree;
    -- We mainly remember F.cache.TorsionFree (as a quotient module of M) for use in SheafMaps.m2.
    F.cache.GlobalSectionLimit = max(0, p);
    -- this is the module Gamma_* M, as a module over the original ring R.
    if p <= 0 then (
	F.cache.SaturationMap = quot;
	-- The map from M to the simplified R-module N.
	-- Since we keep the same module N in this case, Macaulay2 remembers
	-- any earlier calculations done for N, such as Ext calculations.
	F.cache.SaturatedLift = N'; -- This is N as an S-module.
	return
	if complete
	then (-infinity, -infinity, N)
	else (p + bound, p + bound, N));
    G := minimalPresentation target(
	-- an ideal in S generated by variables x_i ^ ceiling(p / d_i).
	-- TODO: substitute with appropriate irrelevant ideal here
	-- TODO: rename this
	Bp := ideal apply(gens S, degrees S, (x, d) -> x ^ (ceiling(p / d#0)));
	-- This is the R-linear map Ip **_S R -> R.
	-- So R-linear maps from the first module to N (which may be M/M_tors)
	-- are equivalent to S-linear maps from Ip to N' = (N as an S-module).
	-- TODO: this inducedMap is faster over S than R, but perhaps
	-- the inclusion of the ideal should be a cached canonicalMap
	inc := inducedMap(module S, module Bp) ** R;
	iso := inducedMap(Hom(R, N), N);
	-- we compute the map N -> Gamma_* F as a limit by
	-- applying Hom(-, N) to the sequence above
	-- cf. the limit from minimalPresentation hook
	-- and emsbound in NormalToricVarieties/Sheaves.m2
	phi := Hom(inc, N, MinimalGenerators => true) * iso);
    -- Note that phi: N -> G is always injective, because we have arranged for N to be m-torsion-free.
    -- If phi is surjective, hence an isomorphism, then we will take the output to be the original module N rather than G,
    -- so that Macaulay2 remembers any earlier calculations done for N, such as Ext calculations.
    if isSurjective phi then (
	F.cache.SaturationMap = quot; -- The map from M to the simplified R-module N.
	F.cache.SaturatedLift = N';
	return
	if complete
	then (-infinity, -infinity, N)
	else (bound, bound, N));
    -- now we compute the center map in the sequence, where m is the maximal ideal of S:
    -- 0 -> HH^0_m(M) -> M -> Gamma_* M -> HH^1_m(M) -> 0
    iota := inverse G.cache.pruningMap; -- map from Gamma_* M to its minimal presentation
    -- quot is the map from M to N (which may be M/M_tors).
    F.cache.SaturationMap = if p <= 0 then iota * quot else iota * phi * quot;
    F.cache.SaturatedLift = minimalPresentation liftModule G;
    if complete then (-infinity, -infinity, G)
    else (bound, bound, G))

-----------------------------------------------------------------------------
-- cohomology
-----------------------------------------------------------------------------
-- TODO: add hooks for X not finite type over k?

-- HH^p(X, OO_X)
cohomology(ZZ,          SheafOfRings) := Module => opts -> (p,    O) -> cohomology(p, variety O, O^1, opts)
cohomology(ZZ, Variety, SheafOfRings) := Module => opts -> (p, X, O) -> cohomology(p,         X, O^1, opts)

-- HH^p(X, F(>=b))
--
cohomology(ZZ,                    SumOfTwists) := Module => opts -> (p,    S) -> cohomology(p, variety S, S, opts)
cohomology(ZZ, ProjectiveVariety, SumOfTwists) := Module => opts -> (p, X, S) -> (
    -- HH^p(X, F(>=b)) computes the cohomology of a coherent sheaf F on a closed subspace X = Proj R
    -- of a weighted projective stack as a graded R-module that maps to HH^p(X, F(*)) in degrees >= b.
    -- The map will be an isomorphism in degrees at least b, but if HH^p(X, F(*)) is bounded below,
    -- then HH^p(X, F(*)) is returned. If p = 0, the map is stored in F.cache.SaturationMap.
    --
    -- H.cache.interval stores the confidence degrees (b0, b1), which indicates that
    -- H maps to HH^p(X, F(*)) isomorphically in degrees >= b0 and surjectively in degrees >= b1.
    -- e.g. b0 == -infinity if and only if the output is HH^p(X,F(*)) and is bounded below.
    --
    -- The distinction between a X as a stack and its associated coarse moduli space f: X -> V
    -- does not matter for computing cohomology. Indeed, a finitely generated graded module M
    -- determines a coherent sheaf M~ on X, and it is a fact that HH^i(X, M~) = HH^i(V, f_*(M~)).
    -- Twists are interpreted by tensoring with the line bundles O(c) on the stack,
    -- which ensures that HH^i(X, M~(c)) = HH^i(X, M(c)~).
    --
    checkVariety(X, S);
    strategy := opts.Strategy;
    verboseLog := if debugLevel > 0 then printerr else identity;
    --
    (F, bound) := (S#0, S#1#0); -- S = F(>=bound)
    -- check if F is a twist of another coherent sheaf
    -- and concentrate cached info in the base twist
    if F.cache.?BaseTwist then (
        verboseLog "computing cohomology of the base twist";
        (F0, twist) := F.cache.BaseTwist; -- F = F0(twist)
        H0 := cohomology(p, X, F0(>= bound + twist#0), opts);
        H  := H(twist);
        if p == 0 then
        F.cache.GlobalSectionLimit  = F0.cache.GlobalSectionLimit;
        F.cache.TorsionFree         = F0.cache.TorsionFree(twist);
        F.cache.SaturationMap       = F0.cache.SaturationMap(twist);
        F.cache.SaturatedLift       = F0.cache.SaturatedLift(twist);
        H.cache.Degrees = (b0, b1) := H0.cache.Degrees - (twist, twist);
        verboseLog("computed cohomology module maps to the cohomology ",
            "isomorphically in degrees >= ", b0, " and surjectively in degrees >= ", b1, ".");
        return H);
    --
    -- handle caching and strategies
    F.cache.HH ??= new MutableHashTable;
    -- TODO: HH^2(F(>=-3)) is in principle determined by HH^2(F(>=-10)),
    -- but it would give a simpler module, which might be preferable.
    -- should we recompute cohomology in this case?
    if F.cache.HH#?p then (
        H = F.cache.HH#p;
        (b0, b1) = H.cache.interval ?? ({bound}, {bound});
        if {bound} < b0 then (
            verboseLog "previously cached cohomology is insufficient";
            remove(F.cache.HH, p)));
    --
    H = F.cache.HH#p ??= runHooks((HH, ZZ, ProjectiveVariety, SumOfTwists),
        (opts, p, X, F, bound), Strategy => strategy);
    if H === null then if strategy === null
    then error("no applicable strategy found for given sheaf cohomology computation")
    else error("assumptions for sheaf cohomology strategy ", toString strategy, " are not met");
    --
    -- check if the lower limit of the confidence interval is good enough
    (b0, b1) = H.cache.interval ?? ({bound}, {bound});
    verboseLog("computed cohomology module maps to the cohomology ",
        "isomorphically in degrees >= ", b0, " and surjectively in degrees >= ", b1, ".");
    H)

addHook((HH, ZZ, ProjectiveVariety, SumOfTwists), Strategy => Ext,
    -- strategy HH^p(X, F(>=b), Strategy => Ext) uses extensions.
    -- that is, HH^p(X, F(a)) = HH^p(X, Hom_X(OO_X(-a), F)) = Ext^p(OO_X(-a), F)
    -- for p = 0 this is similar to the default strategy, but
    -- for p > 0 there can be big differences in speed.
    -- TODO: add hypotheses
    (opts, p, X, F, bound) -> cohomologyByExtensions(p, F, bound))

addHook((HH, ZZ, ProjectiveVariety, SumOfTwists), Strategy => Local,
    -- strategy HH^p(F(>=b), Strategy => Local) uses local duality.
    (opts, p, X, F, bound) -> if 0 < p then cohomologyByLocalDuality(p+1, F, bound))

addHook((HH, ZZ, ProjectiveVariety, SumOfTwists), Strategy => Default,
    -- Note: other strategies must be called directly
    (opts, p, X, F, bound) -> (
        F.cache.HH#p ??= if p == 0
        then twistedGlobalSectionsModule(F, bound)
        else HH^(p+1)(module F,   Degree => bound))
    )

-----------------------------------------------------------------------------

protect Direct
protect DirectNonPrint
protect NonPrint

-- HH^p(F(>=b), Strategy => Ext)
--
cohomologyByExtensions = (opts, p, X, F, bound) -> (
    -- when p > 0 this strategy uses a completely different algorithm based on
    -- HH^p(X, M^~(a)) = Ext^p_R(res B, M)_a where B = (x_0^b_0,...,x_(n-1)^b_(n-1)).
    -- This can be much faster or much slower than the default algorithm.
    --
    -- Note that while the default algorithm gives the complete cohomology module
    -- when the cohomology is bounded below, this algorithm only returns a module
    -- which is correct in degrees >= b, and may take more time if you ask for a
    -- larger range of twists, even if the cohomology is zero in part of that range.
    --
    -- One choice for speed is using the complete intersection ideal B defined above,
    -- rather than the truncations R_(>=j) which have more complicated free resolutions.
    -- This function also uses a bound based on the minimal free resolution of M,
    -- which is sharper than just using the regularity of M.
    --
    -- It might be reasonable to first run hh^i(X, F(*))), which computes the Hilbert series
    -- the cohomology module in all degrees, and then (if you want cohomology as a module)
    -- apply this function starting in the lowest degree where the cohomology is not zero.
    --
    M := module F;
    R := assertWeightedZZGraded ring M;
    quot := torsionFreeMap F; -- the map M -> N = M/M_tors
    N0 := torsionFreeLift F;  -- kill m-torsion and pushforward to S
    N := target quot;
    S := ring N0;
    degs := degrees S; -- This is a list of the form {{1},{9},{15},{22}}, say.
    n := numgens S; -- 1 + dim Proj S
    w := sum degs;
    -- TODO: do we need more than poincare N0?
    bettitable := betti freeResolution(N0, LengthLimit => n - p);
    -- The Betti numbers of the minimal free resolution
    -- ... -> F_1 -> F_0 -> N0 -> 0, correct up to F_{n-p}.
    nonzerokeys := keys bettitable;
    relevantpart0 := select(nonzerokeys, (i,d,h) -> (i == n-p));
    relevantpart1 := select(nonzerokeys, (i,d,h) -> (i == n-p-1));
    maxdeg0 := max apply(relevantpart0, (i,d,h) -> h); -- The maximum degree of a generator of F_{n-p}.
    maxdeg1 := max apply(relevantpart1, (i,d,h) -> h); -- The maximum degree of a generator of F_{n-p-1}.
    maxdeg := max(maxdeg0, maxdeg1); -- The maximum degree of a generator of F_{n-p} or F_{n-p-1}.
    -- It is -infinity if those two free modules are zero.
    -- Then the map Ext^c_S(Bp,N0)_a -> HH^c(X, F(a)) is an isomorphism for all a >= bound, for any homogeneous ideal Bp in S that lives
    -- in degrees > maxdeg - w - bound and that contains a power of the irrelevant ideal. So, define j by:
    -- maxdeg+1-w-bound;
    -- Previous versions took j to be reg_{Symonds}(N0) + 1 - p - bound = reg_{Macaulay2}(N0) - w + n + 1 - p - bound,
    -- which gave a weaker result.
    -- Also, for this choice of j and hence Bp, the map from Ext^c to HH^c is surjective
    -- in degrees a >= bound - (maxdeg-maxdeg1). That adds information if maxdeg1 < maxdeg0.
    -- Finally, if p = 0, then (because N is m-torsion-free) the map from Ext^c to HH^c is injective in all degrees.
    -- So, in this case, we can instead take j = maxdeg1+1-w-bound,
    -- and then the map is surjective (hence an isomorphism) in degrees >= bound.
    j := if p == 0
    then maxdeg1 + 1 - w - bound
    else maxdeg  + 1 - w - bound;
    --
    output := null;
    -- We need a homogeneous S-ideal Bp whic is concentrated in degrees at least j
    -- and that contains a power of the irrelevant ideal (x_0,...,x_(n-1)).
    -- The following seems like the best choice for efficiency:
    H := if j <= 0
    -- take Bp = S, which gives the following results
    then if p == 0 then N else R^0
    -- take Bp = ideal( x_i ^ ceiling(j / d_i) )
    else (
	Bp := ideal apply(gens S, degs,
	    (x, d) -> x ^ (ceiling(j / d#0)));
	--
	-- TODO: can we combine these?
	if p == 0 then (
	    -- If p = 0, then we know that N < output < HH^0(X, F(*)).
	    -- Thus the output is always at least as good an approximation to HH^0(X, F(*)) as N was,
	    -- and so we will cache it (unless it is equal to N).
	    inc := inducedMap(module S, module Bp) ** R;
	    -- This is the R-linear map Bp **_S R -> R. (So R-linear maps from the first module
	    -- to N are equivalent to S-linear maps from Bp to N0 (which is N as an S-module).)
	    iso := inducedMap(Hom(R, N), N);
	    -- we compute the map N -> Gamma_* F as a limit by
	    -- applying Hom(-, N) to the sequence above
	    phi := Hom(inc, N, MinimalGenerators => true) * iso;
	    -- Thus phi: N -> G = Hom_S(Bp, N0), viewed as an R-module.
	    -- Note that phi: N -> G is always injective, because we have arranged for N to be m-torsion-free.
	    -- If phi is surjective, hence an isomorphism, then we will take the output to be the original module N rather than G,
	    -- so that Macaulay2 remembers any earlier calculations done for N, such as Ext calculations.
	    if isSurjective phi
	    then output = N
	    else (
		output = target phi; -- This is Hom_S(Bp, N0), viewed as an R-module. It may be given as a subquotient, which seems fine.
		iota := inverse output.cache.pruningMap; -- the map from G to its minimal presentation, output.
		-- quot is the map from M to N (which may be M/M_tors).
		F.cache.SaturationMap = iota * phi * quot;
		F.cache.SaturatedLift = minimalPresentation liftModule output)
	    )
	else (
	    -- Following the advice of the Macaulay2 documentation, for c = p > 0,
	    -- we compute the graded module Ext^{c+1}(S/Bp,N0), rather than Ext^c(Bp,N0) (which is isomorphic).
	    --
	    -- Also, since we take Bp to be the complete intersection ideal Bp=(x_0^b_0,...,x_(n-1)^b_(n-1)), the minimal free resolution
	    -- of S/Bp is a Koszul complex, which is self-dual; so Ext^{c+1}_{S}(S/Bp, N0) = Tor_{n-c-1}^{S}(S/Bp, N0)(sum |x_i^b_i|),
	    -- keeping track of the grading. That avoids dualizing and hence should be faster.
	    -- Furthermore, since we compute Tor by resolving S/Bp, it is manifestly computed by a complex of R-modules
	    -- (using N in place of N0). That avoids having to tensor with R at the end of the calculation, which can be slow,
	    -- because Macaulay2 would have to convert the Ext module from a subquotient module to a quotient.
	    -- Koszul complex on the row matrix (x_0^b_0,...,x_(n-1)^b_(n-1)) over R.
	    K := koszulComplex(gens Bp, Concentration => (n-p-2, n-p));
	    Torshift := first sum degrees Bp; -- This is sum |x_i^b_i|.
	    output = R^{Torshift} ** HH_(n-p-1)(K ** N);
	    );
	output
	);
    -- For c = p > 0, the map from the output module to HH^c(X,F(*)) is an isomorphism in degrees >= bound+min(j,0)
    -- and surjective in degrees >= bound-(maxdeg-maxdeg1)+min(j,0). Moreover, when p > 0 and j <= 0,
    -- the output module is 0, and so surjectivity implies isomorphism in that case.
    -- For p = 0, the map from the output module to HH^0(X,F(*)) is always injective, and it is surjective
    -- (hence an isomorphism) in degrees >= bound+min(j,0).
    -- We draw the following conclusions.
    --
    -- output = truncate(bound, output);
    -- One might prefer to truncate at bound. As it is, we output a module that is only known
    -- to be correct in degrees at least bound. But this seems OK, especially since we print an explanation.
    -- Truncating can be slow.
    --
    -- the confidence degrees (b0, b1) indicates that H maps to HH^p(X, F(*))
    -- isomorphically in degrees >= b0 and surjectively in degrees >= b1.
    H.cache.interval = (
	if j === -infinity or (p > 0 and j <= 0 and maxdeg1 === -infinity)
	then ({-infinity}, {-infinity})
	else (
	    if maxdeg1 >= maxdeg0 or p == 0
	    then ({bound + min(j,0)}, {bound + min(j,0)})
	    else ( -- Now we in particular have p > 0.
		if maxdeg1 === -infinity
		then ({bound + min(j,0)}, {-infinity})
		else ({bound + min(j,0)}, {bound - (maxdeg - maxdeg1) + min(j,0)}))));
    H)

-- HH^p(X, F)
--
-- As well as caching individual cohomology groups, we try to fix the modules involved,
-- so that Macaulay2 automatically remembers Ext calculations done with them.
cohomology(ZZ,                    CoherentSheaf) := Module => opts -> (p,    F) -> cohomology(p, variety F, F, opts)
cohomology(ZZ,     AffineVariety, CoherentSheaf) := Module => opts -> (p, X, F) -> (
    checkVariety(X, F);
    if p == 0 then module F else (ring X)^0)
cohomology(ZZ, ProjectiveVariety, CoherentSheaf) := Module => opts -> (p, X, F) -> (
    -*
    checkVariety(X, F);
    F.cache.HH   ??= new MutableHashTable;
    if F.cache.HH#?p then return F.cache.HH#p;
    -- TODO: only need basis(0, G) in the end, is this too much computation?
    G := if p == 0 then twistedGlobalSectionsModule(F, 0) -- HH^0 F(>=0)
    else (
	-- pushforward F to PP^n
	M := liftModule module F;
	S := ring M;
	-- TODO: both n and w need to be adjusted for the multigraded case
	n := dim S-1;
	w := S^{-n-1};
	-- using Serre duality for coherent sheaves on schemes with mild
	-- singularities, Cohen–Macaulay schemes, not just smooth schemes.
	-- TODO: check that X is proper (or at least finite type)
	Ext^(n-p)(M, w));
    k := coefficientRing ring X;
    F.cache.HH#p = k^(rank source basis(0, G)))
    *-
    checkVariety(X, F);
    F.cache.HH   ??= new MutableHashTable;
    if F.cache.HH#?p then return F.cache.HH#p;
    k := coefficientRing ring variety F;
    F.cache.HH#p = k^(hh^p(F)))

-- HH^p(OO_X)
cohomology(ZZ,          SheafOfRings) := Module => opts -> (p,    O) -> cohomology(p, O^1, opts)
cohomology(ZZ, Variety, SheafOfRings) := Module => opts -> (p, X, O) -> cohomology(p, X, O^1, opts)

-- HH^p(F(>=b), Strategy => Local)
--
-- This is essentially the local cohomology algorithm from m2/local.m2,
-- but adjusted to interpret the truncation option correctly in the weighted case
-- TODO: remove m2/local.m2 and use this for HH^ZZ(Module), and add HH_B^ZZ(Module)
cohomologyByLocalDuality = (opts, X, i, F, bound) -> (
    -- This computes the local cohomology HH_m^i(M) at the maximal ideal m.
    -- When i > 0, this is a graded R-module which maps to HH^i(X, F(*))
    -- and is an isomorphism in degrees >= b0 (with b0 <= the input b)
    -- and is a surjection in degrees >= b1.
    if bound == -infinity then notImplemented();
    R := ring module F;
    M := torsionFreeLift F;
    S := ring M;
    n := numgens S; -- 1 + dim Proj S
    -- Macaulay2's function for local cohomology (corrected below)
    -- gives accurate answers even if R has generators in different degrees;
    -- but the degree from which it starts computing is off by sigma = w - n,
    -- which explains the change made here to the definition of the dualizing module w.
    -- If HH^i(X, F(*)) is bounded below, however, then that is irrelevant;
    -- the output is the whole cohomology, regardless of the given bound.
    w := S.cache.Dualizing ??= S^{-sum degrees S};
    -- We fix the dualizing module w, as a graded R-module.
    -- As a result, Macaulay2 automatically remembers Ext^(n-i)(M, w),
    -- in case another function has computed that module earlier.
    E := minimalPresentation Ext^(n-i)(M, w);
    --
    -- initial confidence interval
    bounds := ({bound}, {bound});
    result := if dim E <= 0 then (
	bounds = ({-infinity}, {-infinity});
	Ext^n(E, w))
    else (
	topdegree := - first min degrees E;
	-- In this case, HH^i(X, F(*)) is unbounded below,
	-- and topdegree is the top degree in which it is not zero.
	if bound > topdegree then return (
	    output := R^0;
	    output.cache.interval = ({topdegree + 1}, {topdegree + 1});
	    output);
	truncatedDual(E, bound));
    --
    output = result ** R; -- This can be somewhat slow, because computing the tensor product requires Macaulay2
    -- to find a presentation of the S-module "result" (which will usually be given as a subquotient).
    -- This happens even though "result" is already the S-module underlying an R-module.
    -- At least Macaulay2's tensor product algorithm avoids finding a presentation in the case when S = R.
    output.cache.interval = bounds;
    output)

-- Given a graded module E over a (possibly nonstandard) graded polynomial ring S over kk,
-- return a graded module that agrees with the kk-dual of E in degrees >= -e.
-- cf. truncatedDual in m2/local.m2.
-- The output of this version may be nonzero in degrees less than e for speed.
truncatedDual = (E, e) -> (
    -- find kk-dual of E, truncated in degrees >= e.
    -- depends on truncate methods
    S := ring E;
    n := numgens S;
    degs := degrees S;
    -- 
    -- We will define E1 to be a quotient module of E that has finite length
    -- and that agrees with E in degrees at most -e.
    degsE := degrees E; -- The degrees of the generators of E, in the form {2,3,...}.
    gensToKill := positions(degsE, i -> i  > {-e}); -- We will kill these generators in E1.
    gensToKeep := positions(degsE, i -> i <= {-e});
    rels1 := apply(gensToKill, i -> E_i);
    -- multiply each generator of E with degree d_i <= -e
    -- by each variable x_j to the power roundup((-e+1-d_i)/a_j).
    rels2 := flatten apply(gensToKeep,
	i -> apply(n, j -> ((S_j)^-((degsE#i#0 - 1 + e) // degs#j#0) * E_i)));
    -- At the moment, just writing E/(rels1|rels2) would typically not
    -- be homogeneous, hence this more complex formulation.
    E1 := minimalPresentation (E / image map(E, , matrix(rels1 | rels2)));
    --
    w := S.cache.Dualizing ??= S^{-sum degs};
    Ext^n(E1, w))

-----------------------------------------------------------------------------
-- Module of twisted global sections Γ_*(F)
-----------------------------------------------------------------------------

-- This is an approximation of Gamma_* F, at least with an inclusion from Gamma_>=0 F
-- TODO: optimize caching: if HH^0(F(>=b)) is cached above, does this need to be cached?
minimalPresentation SheafOfRings  := prune SheafOfRings  := SheafOfRings  => opts -> identity
minimalPresentation CoherentSheaf := prune CoherentSheaf := CoherentSheaf => opts -> (
    F -> F.cache#(symbol minimalPresentation => opts) ??= tryHooks(
        (minimalPresentation, CoherentSheaf), (opts, F),
        (opts, F) -> sheaf prune(module F, opts)))

addHook((minimalPresentation, CoherentSheaf), Strategy => ProjectiveVariety,
    (opts, F) -> if isProjective variety F then (
        -- HH^0(F(>=0)) uses twistedGlobalSectionsModule and caches SaturationMap.
        -- We want the target of the saturation map gamma: M -> Gamma_(d >= 0)(X, F(d)),
        -- so we first check if SaturationMap is cached, and if not we call HH^0(F(>=0))
        if not F.cache.?SaturationMap then HH^0(F(>=0), Degree => NonPrint);
        gamma := sheaf(F.variety, F.cache.SaturationMap);
        G := target gamma; -- G = HH^0(F(>=0))
        G.cache.pruningMap = gamma;
        G)
    )

addHook((minimalPresentation, CoherentSheaf), Strategy => "SerreTwist",
    (opts, F) -> if isProjective variety F and F.cache.?BaseTwist then (
        (F0, twist) := F.cache.BaseTwist; -- F = F0(twist)
        G := (G0 := minimalPresentation(F0, opts))(twist);
        -- these would normally be set by twistedGlobalSectionsModule,
        -- so we transfer them in the twisted sheaf for use in SheafMaps.m2
        F.cache.GlobalSectionLimit = F0.cache.GlobalSectionLimit;
        F.cache.TorsionFree        = F0.cache.TorsionFree(twist);
        G.cache.pruningMap         = G0.cache.pruningMap(twist);
        G)
    )

-----------------------------------------------------------------------------
-- Projective bundles
-----------------------------------------------------------------------------
-- TODO: add isVectorSpace, then given a vector space V with basis elements
-- V_1 .. V_n support defining PP(V) = Proj Sym V = Proj kk[V_1..V_n].

-- TODO: is this correct?
symmetricAlgebra CoherentSheaf := Ring => opts -> F -> symmetricAlgebra(HH^0 F(>=0), opts)
-- TODO: is the dual right?
-- TODO: add isLocallyFree and make sure F is locally free first?
ProjectiveSpace(CoherentSheaf) := ProjectiveVariety => F -> tryHooks((ProjectiveSpace, CoherentSheaf), F,
    F -> Proj flattenRing(symmetricAlgebra dual F, Result => Thing))

-----------------------------------------------------------------------------
-- Sheaf Hom and Ext
-----------------------------------------------------------------------------

sheafHom = method(TypicalValue => CoherentSheaf, Options => options Hom)
sheafHom(SheafOfRings, SheafOfRings)  :=
sheafHom(SheafOfRings, CoherentSheaf) :=
sheafHom(CoherentSheaf, SheafOfRings)  :=
sheafHom(CoherentSheaf, CoherentSheaf) := CoherentSheaf => opts -> (F, G) -> (
    assertSameVariety(F, G); sheaf(variety F, Hom(module F, module G, opts)))

sheafExt = new ScriptedFunctor from {
    superscript => i -> new ScriptedFunctor from {
	-- sheafExt^1(F, G)
	argument => X -> applyMethod''(sheafExt, functorArgs(i, X))
	},
    argument => X -> applyMethod''(sheafExt, X)
    }

sheafExt(ZZ, SheafOfRings, SheafOfRings)  :=
sheafExt(ZZ, SheafOfRings, CoherentSheaf) :=
sheafExt(ZZ, CoherentSheaf, SheafOfRings)  :=
sheafExt(ZZ, CoherentSheaf, CoherentSheaf) := CoherentSheaf => options Ext.argument >> opts -> (i, F, G) -> (
    assertSameVariety(F, G); sheaf(variety F, Ext^i(module F, module G, opts)))

-----------------------------------------------------------------------------
-- Global Ext via efficient truncation.
-- We build on ideas from the algorithm in:
-- Gregory G. Smith, Computing global extension modules,
-- Journal of Symbolic Computation, 29 (2000) 729-746.
-- See documentation of Ext^ZZ(CoherentSheaf,CoherentSheaf) for examples.
-----------------------------------------------------------------------------

-- Given an integer p, a complex D = Hom(res M, N), where M and N are torsion-free
-- graded modules over a graded ring R = S/I representing sheaves on X = Proj R,
-- and a bound b, return an integer e (or -infinity) such that:
--   if B is any ideal in the polynomial ring S such that B contains S in
--   sufficiently high degrees and B is in graded degree >= e, then the map
--    Ext^p_S(B, D)_(>=b) -> HH^p(X, D~(>=b)) is an isomorphism.
--
-- The function returns a sequence (b0, b1, e), meaning that the output
-- of Ext^p_X(F, D^~(>=b)) (using the number e) will be a module that maps isomorphically
-- to Ext^p_X(F, D^~(*)) in degrees >= b0 and surjectively in degrees >= b1.
degreeBound = (p, D, bound) -> (
    D' := liftComplex D; -- pushforward of the complex to the polynomial ring S
    S := assertWeightedZZGraded ring D';
    n := numgens S; -- 1 + dim Proj S
    w := first sum degrees S;
    -- TODO: can we use the poincare polynomial of the terms somehow?
    bettitable := betti freeResolution(D', LengthLimit => n - p - min D');
    -- The Betti numbers of the minimal free resolution
    -- ... -> F_(1 + min D') -> F_(min D') -> 0 of D' over S, correct up to F_{n-p}.
    nonzerokeys := keys bettitable;
    relevantpart0 := select(nonzerokeys, (i,d,h) -> (i == n-p));
    relevantpart1 := select(nonzerokeys, (i,d,h) -> (i == n-p-1));
    maxdeg0 := max apply(relevantpart0, (i,d,h) -> h); -- The maximum degree of a generator of F_{n-p}.
    maxdeg1 := max apply(relevantpart1, (i,d,h) -> h); -- The maximum degree of a generator of F_{n-p-1}.
    maxdeg := max(maxdeg0, maxdeg1); -- The maximum degree of a generator of F_{n-p} or F_{n-p-1}.
    -- It is -infinity if those two free modules are zero.
    -- Then the map Ext^c_S(B, D')_a -> HH^c(X, D~(a))
    -- is an isomorphism for all a >= bound, for any homogeneous ideal B in S that lives
    -- in degrees > maxdeg - w - bound and that contains a power of the irrelevant ideal.
    -- So we define e = maxdeg + 1 - w - bound.
    -- Also, for this choice of e and hence Bp, the map from Ext^c to HH^c is surjective
    -- in degrees >= bound - (maxdeg - maxdeg1). That adds information if maxdeg1 < maxdeg0.
    --
    -- Previous versions took e to be
    --   reg_{Symonds}(N0) + 1 - p - bound = reg_{Macaulay2}(N0) - w + n + 1 - p - bound,
    -- which gave a weaker result.
    --
    -- TODO: does this apply also to Hom complexes?
    -- if p = 0, then (because N is m-torsion-free) the map
    -- from Ext^0 to HH^0 is injective in all degrees.
    -- in this case, we can instead take e = maxdeg1 + 1 - w - bound,
    -- and then the map is surjective (hence an isomorphism) in degrees >= bound.
    e := if p == 0
    then maxdeg1 + 1 - w - bound
    else maxdeg  + 1 - w - bound;
    -- Note: if e = -infinity, both of these become -infinity
    b0 := bound + min(e,0);
    b1 := bound + min(e,0) - (maxdeg - maxdeg1);
    (b0, b1, e))

    -- if e === -infinity or (p > 0 and e <= 0 and maxdeg1 === -infinity)
    -- then ({-infinity}, {-infinity})
    -- else (
    -- 	if maxdeg1 >= maxdeg0 or p == 0
    -- 	then ({bound + min(e,0)}, {bound + min(e,0)})
    -- 	else (
    --      -- Now we in particular have p > 0.
    -- 	    if maxdeg1 === -infinity
    -- 	    then ({bound + min(e,0)}, {-infinity})
    -- 	    else ({bound + min(e,0)}, {bound + min(e,0) - (maxdeg - maxdeg1)})))


protect TruncateDegree
-- Ext^m(F, G(>=b)), for coherent sheaves F and G over a projective scheme
-- (or weighted projective stack) X = Proj R2 and an integer m,
-- returns a graded R2-module that agrees with Ext^m_X(F, G(*)) in degrees >= b.
-- If we introduce commands for Ext between complexes of sheaves, this function may need to watch
-- for the case where F is a sheaf and G is a complex of sheaves. The option "MinimalGenerators => NonPrint"
-- avoids printing; instead, it returns a sequence (b0, b1, M0) with M0 a module
-- that maps to Ext^m_X(F, G(*)), isomorphically in degrees >= b0 and surjectively in degrees >= b1.
Ext(ZZ, SheafOfRings,  SumOfTwists) := Module => opts -> (m, O, S) -> Ext^m(O^1, S, opts)
Ext(ZZ, CoherentSheaf, SumOfTwists) := Module => opts -> (m, F, S) -> (
    (G, bound) := (S#0, S#1#0); -- Here G should be a coherent sheaf.
    -- TODO: turn this into a hook prehaps?
    if F.cache.?BaseTwist or G.cache.?BaseTwist then ( -- Here F or G was defined as a twist of another sheaf,
	-- say F = F0(c) and G = G0(d). We reduce to the calculation for F0 and G0, to take advantage of caching.
	(F0, twistF) := if F.cache.?BaseTwist then F.cache.BaseTwist else (F, {0});
	(G0, twistG) := if G.cache.?BaseTwist then G.cache.BaseTwist else (G, {0});
	E0 := Ext^m(F0, G0(>= bound + (twist := twistG - twistF)), opts);
	E := E0(twist);
	E.cache.interval = E0.cache.interval - ({twist}, {twist});
	-- TODO: double check that this doesn't need to be twisted
	E.cache.TruncateDegree = E0.cache.TruncateDegree)
    else (
	-- Now the sheaves F and G were not defined as twists.
	N := target torsionFreeMap G; -- A simplified module that represents G.
	R2 := assertWeightedZZGraded ring N;
	degs := degrees R2;
	n := #degs;
	M := target torsionFreeMap F;
	Mres := freeResolution(M, LengthLimit => m+1);
	-- TODO: concentration for Hom?
	HomMN := Hom(Mres, N);
	(b0, b1, e) := degreeBound(m, HomMN, bound);
	-- We need to construct a homogeneous ideal Bp in the graded polynomial ring S that is concentrated
	-- in degrees at least e and that contains a power of the irrelevant ideal (x_0,...,x_(n-1)). The following
	-- seems like the best choice for efficiency. We tensor the Koszul complex that resolves this ideal
	-- over S with R2, so we can work entirely with R2-modules.
	-- In this case, we can take Bp = S, which gives the following.
	E = if e <= 0 then E = HH^m(HomMN)
	else (
	    -- Here e > 0.
	    M2gens := apply(n, i -> ((R2_i)^-((-e)//degs#i#0)));
	    -- That is, the list M2gens consists of each variable x_i to the power b_i := roundup(e/a_i).
	    koszulRes := (koszulComplex(M2gens, Concentration => (max(0, n-m-2), n-1)))[n-1];
	    -- Thus koszulRes is in homological degrees 0, -1, ..., -min(n-1, m+1).
	    Torshift := first sum(M2gens, degree); -- This is sum_i a_i b_i, which is at least ne.
	    (HH^m(koszulRes ** HomMN))(Torshift));
	-- The map Ext^m_R2(M2 tensor_(S) R2, Hom(Mres, N)) -> Ext^m_X(F, G(*)) is an isomorphism in graded degrees >= bound,
	-- and so we will return the first module. We don't truncate or prune it; truncating can be slow,
	-- and this module will often be simpler (in terms of generators and relations) than a truncation
	-- of it. The user is told what the output means. Eventually, we could add information
	-- about surjectivity of the map in some degrees, as in cohomologyByExtensions.
	E.cache.TruncateDegree = e);
    -- Here e is the degree at which the polynomial ring S was truncated to form M2;
    -- this may be used later for computing the Yoneda extension.
    --
    E)

Ext(ZZ, SheafOfRings, SheafOfRings)  :=
Ext(ZZ, SheafOfRings, CoherentSheaf) :=
Ext(ZZ, CoherentSheaf, SheafOfRings)  :=
Ext(ZZ, CoherentSheaf, CoherentSheaf) := Module => opts -> (n, F, G) -> (
    E := (Ext^n(F, G(>=0), opts ++ { MinimalGenerators => NonPrint }))#2;
    -- With the NonPrint option, Ext returns a sequence (b0,b1,M0), and we just want the module M0.
    k := coefficientRing ring E;
    V := k^(hilbertFunction(0, E));
    V.cache.formation = FunctionApplication { Ext, (n, F, G) };
    V.cache.TruncateDegree = E.cache.TruncateDegree;
    V)

-----------------------------------------------------------------------------
-- hh: Hodge decomposition
-----------------------------------------------------------------------------
-- TODO: HodgeTally for pretty printing the Hodge diamond

hhOptions = new OptionTable from {
    Degree => 0
    }

-- Modeled on Ext, defined as a ScriptedFunctor in complexes.m2.
hh = new ScriptedFunctor from {
    superscript => i -> new ScriptedFunctor from {
	-- hh^i(F), hh^i(F, a, b), and so on
	argument => hhOptions >> opts -> X -> applyMethodWithOpts''(hh, functorArgs(i, X), opts)
	},
    argument => hhOptions >> opts -> X -> applyMethodWithOpts''(hh, X, opts)
    }

-- using Hodge symmetry and Serre duality to ease the computation
-- TODO: is minimum necessarily the most efficient?
min'pq = d -> (p,q) -> min{(p,q), (q,p), (d-p,d-q), (d-q,d-p)}

-- The Hodge numbers of a projective variety over a field.
-- By definition, hh^(p,q)(X) = dim HH^q(X, Omega^[p]), in terms of reflexive differentials.
-- If X is smooth, this is the bundle Omega^p of p-forms.
hh(Sequence, ProjectiveVariety) := ZZ => opts -> (pq, X) -> (
    -- cotangentSheaf seems to be the slowest part of this algorithm, so we minimize the exterior powers
    -- (note the definition of hh^(p,q)(X), above).
    (p,q) := if char X == 0 then (min'pq dim X) pq else pq;
    -- Hodge symmetry holds in characteristic 0, but not always in positive characteristic.
    if not X.cache.?hh   then X.cache.hh = new MutableHashTable;
    if X.cache.hh#?(p,q) then X.cache.hh#(p,q) else X.cache.hh#(p,q) = (
	hh^q reflexiveDifferentials(p, X)))

-----------------------------------------------------------------------------

-- This function hh^i(F) computes coherent sheaf cohomology for a coherent sheaf
-- on a closed subspace in a projective space, or more generally in a weighted projective space.
-- If you want to compute cohomology with many twists, the functions hh^i(F,b1,b2) or hh^i(F(*))
-- should be faster than running this program repeatedly.
--
-- This program computes only the dimension of the cohomology, not the cohomology as a vector space.
-- The algorithm uses local duality, as the function HH^i(F(>=a)) does when i>0.
-- It should be faster, in most cases, since it has less to compute.
--
-- The distinction between a WPS as a stack X and its associated coarse moduli space, f: X -> V, does not matter
-- for computing cohomology. Indeed, a finitely generated graded module M determines a coherent sheaf M~
-- on X, and it is a fact that HH^i(X, M~) = HH^i(V, f_*(M~)) for every i. Twists are interpreted by tensoring with
-- the line bundles O(c) on the stack, which ensures that HH^i(X,M~(c)) = HH^i(X,M(c)~).
--
-- We allow the option hh^i(F, Degree => a) to mean hh^i(F(a)). We expect people to use the latter notation,
-- but (in effect) we translate it into the first notation for better caching. (We only need to store information
-- about one sheaf, not its twists.)
--
hh(ZZ, SheafOfRings)  := ZZ => opts -> (cohodeg, O) -> hh^cohodeg(O^1, opts)
hh(ZZ, CoherentSheaf) := ZZ => opts -> (cohodeg, F) -> (
    if not isProjective variety F then ( -- We give a correct answer in the affine case, although that's not our main focus.
	return if cohodeg != 0 then 0
	else if dim F > 0 then infinity
	else degree F);
    if F.cache.?BaseTwist then return hh^cohodeg(F.cache.BaseTwist#0, Degree => opts.Degree + first F.cache.BaseTwist#1);
    -- Thus we reduce to the case where the sheaf F was not defined as a twist.
    R2 := ring module F;
    M := torsionFreeLift F;
    S := assertWeightedZZGraded ring M; -- S is a graded polynomial ring, and M is a simplified S-module that represents the sheaf F.
    -- In particular, we have arranged that M has no m-torsion (where m is the maximal ideal S_(>0)).
    n := numgens S; -- 1 + dim Proj S
    w := S.cache.Dualizing ??= S^{-sum degrees S};
    -- We fix the dualizing module w, as a graded S-module. As a result, Macaulay2 automatically remembers Ext^j(M, w)
    -- (for a number j) in case another function has computed that module earlier.
    output := 0;
    if cohodeg == 0 then (
	output = hilbertFunction(-(opts.Degree), Ext^(n-1-cohodeg)(M, w));
	-- We have arranged that HH^0_m(M) = M_tors is zero, so we don't need to subtract its graded dual Ext^(n-cohodeg)(M, w).
	output = output + hilbertFunction(opts.Degree, M)
	)
    else (
	output = hilbertFunction(-(opts.Degree), Ext^(n-1-cohodeg)(M, w))
	)
    )


-- Replace T by T^{-1} for an element of ZZ[T] (meaning Z[T,T^(-1)]). This replaces
-- "substitute(f, T => T^(-1))", which is slow.
invertvar = (f) -> (
    A := ring f; -- This should be a ring of the form ZZ[T].
    T := A_0; -- This is the variable in that ring.
    coeffseq := (coefficients f)_1; -- This is a column matrix of the coefficients (ignoring the monomials).
    exposeq := flatten exponents f; -- This is the list of exponents of the monomials, in the form {-1, 2, 3}.
    monseq := {apply(exposeq, x -> T^(-x))}; -- This is the list of inverted monomials, in the form {{T^1, T^-2, T^-3}}.
    monomialmatrix := matrix monseq; -- This is the corresponding row matrix.
    (monomialmatrix*coeffseq)_(0,0))


-- This function hh^i(F,b1,b2) computes cohomology with coefficients in a coherent sheaf F
-- on a closed subspace of a weighted projective space, with all twists in an interval [b1,b2].
-- This should be faster than running hh^i(F(b)) separately for many twists b.
-- In some cases, you might prefer hh^i(F(*)), which computes the cohomology in all twists.
--
-- This program computes only the Hilbert series of the cohomology, not the cohomology as a module.
-- The algorithm uses local duality, as in the function HH^i(F(>=b)) that computes the module.
-- It should be faster, in most cases, since it has less to compute.
--
-- The distinction between a WPS as a stack X and its associated coarse moduli space, f: X -> V, does not matter
-- for computing cohomology. Indeed, a finitely generated graded module M determines a coherent sheaf M~
-- on X, and it is a fact that HH^i(X, M~) = HH^i(V, f_*(M~)) for every i. Twists are interpreted by tensoring with
-- the line bundles O(c) on the stack, which ensures that HH^i(X,M~(c)) = HH^i(X,M(c)~).
--
-- We allow the option hh^i(F, b1, b2, Degree => a) to mean hh^i(F(a), b1, b2). We expect people to use the latter notation,
-- but (in effect) we translate it into the first notation for better caching. (We only need to store information
-- about one sheaf, not its twists.)
--
hh(ZZ, SheafOfRings,  ZZ, ZZ) := RingElement => opts -> (cohodeg, O, b1, b2) -> hh^cohodeg(O^1, b1, b2, opts)
hh(ZZ, CoherentSheaf, ZZ, ZZ) := RingElement => opts -> (cohodeg, F, b1, b2) -> (
    checkProjective variety F;
    if F.cache.?BaseTwist then return hh^cohodeg(F.cache.BaseTwist#0, b1, b2, Degree => opts.Degree + first F.cache.BaseTwist#1);
    -- Thus we reduce to the case where the sheaf F was not defined as a twist.
    b1 = opts.Degree + b1; b2 = opts.Degree + b2;
    if instance(F, SheafOfRings) then F = F^1; -- That makes F a CoherentSheaf.
    if not instance(F, CoherentSheaf) or not instance(b1, ZZ) or not instance(b2, ZZ) or b1>b2 then (
	error "the input should be in the form hh^i(F,b1,b2), with F a coherent sheaf and b1 <= b2 integers");
    R2 := ring module F;
    M := torsionFreeLift F;
    S := assertWeightedZZGraded ring M; -- S is a graded polynomial ring, and M is a simplified S-module that represents the sheaf F.
    -- In particular, we have arranged that M has no m-torsion (where m is the maximal ideal S_(>0)).
    A := degreesRing S; -- This is a ring of the form "ZZ[T]" (meaning Z[T,T^(-1)]).
    T := A_0; -- This is the variable in the ring A.
    n := numgens S; -- 1 + dim Proj S
    w := S.cache.Dualizing ??= S^{-sum degrees S};
    -- We fix the dualizing module w, as a graded S-module. As a result, Macaulay2 automatically remembers Ext^j(M, w)
    -- (for a number j) in case another function has computed that module earlier.
    output := 0;
    extmodule1 := Ext^(n-1-cohodeg)(M, w);
    poincare extmodule1; -- This caches the Poincare polynomial of extmodule1 = Ext^(n-1-cohodeg)(M, w),
    -- automatically used by hilbertSeries in what follows.
    if cohodeg == 0 then (
	-- We have arranged that HH^0_m(M) = M_tors is zero, so we don't need to subtract off Ext^n(n-cohodeg)(M, w),
	-- the graded dual to HH^0_m(M).
	output = hilbertSeries(extmodule1, Order => -b1+1);
	output = invertvar(output); -- This is what we want, in degrees at least b1.
	output = output + hilbertSeries(M, Order => b2+1);
	output = part(b1, b2, output) -- We remove terms of degree below b1 or above b2.
	)
    else (
	output = hilbertSeries(extmodule1, Order => -b1+1);
	output = invertvar(output); -- This is the answer, in degrees at least b1.
	output = part(b1, b2, output) -- We remove terms of degree above b2.
	);
    output * T^-(opts.Degree))

-- Compute f1+f2, for f1 a "Divide" (a kind of rational function in a variable T)
-- and f2 a Laurent polynomial in T (as an element of "ZZ[T]" = Z[T,T^(-1)]. The output is again a Divide.
--
add = (f1, f2) -> (
    num1 := value numerator f1;  -- The numerator is a Laurent polynomial in ZZ[T]. (We need to say "value" because
    -- if the numerator is, say, T^3 (rather than a bigger polynomial), then its type may be "Power". To allow different
    -- types of input, we will likewise use value(f2) in what follows.)
    den1 := denominator f1; -- The denominator is a "Product", something like (1-T^2)(1-T^3).
    Divide{num1 + (value f2) * (value den1), den1})

-- Compute f1*f2, for f1 a "Divide" (a kind of rational function in a variable T)
-- and f2 a Laurent polynomial in T (as an element of "ZZ[T]" = Z[T,T^(-1)]. The output is again a Divide.
--
mult = (f1, f2) -> Divide(value(numerator f1)*value(f2),denominator f1)
-- The numerator of f1 is a Laurent polynomial in ZZ[T]. (We need to say "value" because
-- if the numerator is, say, T^3 (rather than a bigger polynomial), then its type may be "Power". To allow different
-- types of input, we likewise use value(f2) here.)
-- The denominator of f1 is a "Product", something like (1-T^2)(1-T^3).

-- Compute the top degree of a series in T and T^{-1}. The input is positiveseries (a "Divide" in a variable T)
-- and negativeseries (a "Divide" in a variable U, viewed as T^(-1)). It is assumed that all terms of negativeseries
-- have degree (in terms of T) below the degree of all terms of positiveseries. The output
-- could be infinity (if positiveseries is not a polynomial in T) or -infinity (if both inputs are 0).
topDegree = (positiveseries, negativeseries) -> (
    posseries := reduceHilbert positiveseries;
    negseries := reduceHilbert negativeseries; -- We simplify the input functions, so we can see whether they are polynomials.
    A := ring numerator negseries;
    U := A_0;
    num := 0;
    if value(denominator(posseries)) != 1 then infinity
    else (  -- Here posseries is a polynomial in T.
	num = numerator(posseries);
	if num != 0 then (degree(num))_0
	else ( -- Here posseries is 0.
	    num = numerator(negseries);
	    if num != 0 then (degree(substitute(num, U => U^(-1))))_0
	    else -infinity -- Here both inputs are zero.
	    )
	)
    )


-- This function (usually called as hh^i(F(*)))
-- computes coherent sheaf cohomology on a closed subspace of a weighted projective space with all twists.
-- For the input hh^i(F(>=b)), the number b is ignored, as the output explains. The base ring should be a field.
-- The output is a sequence listing the infimum of weights c such that HH^i(X,F(c)) is not zero (possibly -infinity or,
-- if the cohomology is zero in all weights, infinity), the supremum of such weights,
-- a rational function in T, and a rational function in U = T^{-1}, such that the sum of these two functions
-- as Laurent series is sum_c h^i(X, F(c)) T^c (the sum over all integers c). (The two functions do not overlap,
-- as formal series in T.)
--
-- This program computes only the Hilbert series of the cohomology, not the cohomology as a module.
-- The algorithm, using local duality, is the same as that used for hh^i(F) and hh^i(F,b1,b2),
-- and slightly different from that used for the module HH^i(F(>=b)). It should be faster, in most cases.
--
-- The distinction between a WPS as a stack X and its associated coarse moduli space, f: X -> V, does not matter
-- for computing cohomology. Indeed, a finitely generated graded module M determines a coherent sheaf M~
-- on X, and it is a fact that HH^i(X, M~) = HH^i(V, f_*(M~)) for every i. Twists are interpreted by tensoring with
-- the line bundles O(c) on the stack, which ensures that HH^i(X,M~(c)) = HH^i(X,M(c)~).
--
-- We allow the option hh^i(F, Degree => a) to mean hh^i(F(a)). We expect people to use the latter notation,
-- but (in effect) we translate it into the first notation for better caching. (We only need to store information
-- about one sheaf, not its twists.)
--
hh(ZZ, SumOfTwists) := Sequence => opts -> (cohodeg, sumoftwists) -> (
    -- Compute HH^{cohodeg}(X,F(a)) for all integers a, as a sum of two Laurent series,
    -- where F is a CoherentSheaf (or a SheafOfRings) on a closed subspace of a weighted projective space.
    F := sumoftwists#0; -- For an input of the form hh^i(F(>=b)), the number b is ignored. Note that,
    -- if F is input as a SheafOfRings, SumOfTwists automatically turns it into a CoherentSheaf; so that's what this function receives.
    if F.cache.?BaseTwist then return hh^cohodeg((F.cache.BaseTwist#0)(*), Degree => opts.Degree + first F.cache.BaseTwist#1);
    -- Thus we reduce to the case where the sheaf F was not defined as a twist.
    R2 := ring module F;
    M := torsionFreeLift F;
    S := assertWeightedZZGraded ring M; -- S is a graded polynomial ring, and M is a simplified S-module that represents the sheaf F.
    -- In particular, we have arranged that M has no m-torsion (where m is the maximal ideal S_(>0)).
    A := degreesRing S; -- This is a ring of the form "ZZ[T]" (meaning Z[T,T^(-1)]).
    T := A_0; -- This is the variable in the ring A.
    U := getSymbol "U"; -- FIXME
    B := newRing(A,Variables=>{U},MonomialOrder=>{MonomialSize=>32,Weights=>{-1},GroupLex=>1,Position=>Up},Inverses=>true);
    U = B_0;
    -- Thus B is a ring of the form "ZZ[U]" (meaning Z[U,U^(-1)]). We think of U as meaning T^(-1).
    n := numgens S; -- 1 + dim Proj S
    w := S.cache.Dualizing ??= S^{-sum degrees S};
    -- We fix the dualizing module w, as a graded S-module. As a result, Macaulay2 automatically remembers Ext^j(M, w)
    -- (for a number j) in case another function has computed that module earlier.
    bottomdeg := 0;
    topdeg := 0;
    positiveseries := 0;
    negativeseries := 0;
    hilb1 := 0; hilb0 := 0;
    extmodule1 := Ext^(n-1-cohodeg)(M, w); -- Macaulay2 automatically caches this, in case Ext was computed earlier or will be computed later.
    poincare extmodule1; -- This caches the Poincare polynomial of extmodule1 = Ext^(n-1-cohodeg)(M, w),
    -- automatically used by hilbertSeries in what follows.
    if cohodeg == 0 then (
	hilbM := mult(hilbertSeries(M, Reduce => true), T^(-opts.Degree)); -- This is a "Divide".
	if dim extmodule1 <= 0 then ( -- Here HH^1_m(M) has finite length, and so (equivalently) HH^0(X,F(*)) is bounded below.
	    -- In this case, we will return the output as a rational function in T.
	    hilb1 = T^(opts.Degree) * value numerator hilbertSeries(extmodule1, Reduce => true); -- The graded dual of this is local cohomology HH^1_m(M),
	    -- which has finite length in the case at hand. So its Hilbert series is a Laurent polynomial.
	    -- We arranged that HH^0_m(M) = M_tors is zero; so we don't need to subtract off the Hilbert series of Ext^n(M, w),
	    -- the graded dual of HH^0_m(M).
	    hilb10 := substitute(hilb1, T => T^(-1)); -- This is a Laurent polynomial in ZZ[T] = Z[T, T^(-1)].
	    positiveseries = add(hilbM, hilb10);
	    negativeseries = 0;
	    if value numerator positiveseries == 0 then (  -- In this case, HH^0(X,M(*)) = 0.
		topdeg = -infinity;
		bottomdeg = infinity)
	    else (
		topdeg = infinity;
		bottomdeg = first min exponents value numerator positiveseries -- This gives the bottom degree of the output, positiveseries.
		)
	    )
	else ( -- Here extmodule1 is not bounded above, and so HH^0(X, F(*)) is not bounded below.
	    -- We will return the output as a rational function in T (of degree >=0 as a series in T) plus a Laurent polynomial in U = T^(-1)
	    -- (of degree >0 as a series in U, hence of degree <0 as a series in T).
	    hilb1 = mult(hilbertSeries(extmodule1, Reduce => true), T^(opts.Degree));
	    pospart1 := hilbertSeries(extmodule1, Order => 1 - opts.Degree) * T^(opts.Degree); -- This is a Laurent polynomial (not series)
	    -- in degrees <= 0, but we'll change T to T^{-1}; hence the name "pospart1".
	    negpart1 := add(hilb1, -pospart1); -- This is a "Divide", the part of hilb1 in T-degree > 0.
	    -- We arranged that HH^0_m(M) = M_tors is zero, so we don't need to subtract off a contribution from Ext^n(M,w),
	    -- the graded dual of HH^0_m(M).
	    pospart10 := substitute(pospart1, T => T^(-1)); -- A polynomial in T.
	    negpart10 := substitute(negpart1, T => U); -- A "Divide" in U, of U-degree >= 1 as a series.
	    negpartM := hilbertSeries(M, Order => opts.Degree) * T^(-opts.Degree); -- The part of hilbM of T-degree < 0.
	    positiveseries = add(hilbM, -negpartM + pospart10);
	    negativeseries = add(negpart10, substitute(negpartM, T => U^(-1)));
	    topdeg = topDegree(positiveseries, negativeseries);
	    bottomdeg = -topDegree(negativeseries, positiveseries) -- This works even if these degrees are infinity or -infinity.
	    )
	)
    else (  -- Here we are computing cohomology in some degree > 0.
	hilb1 = mult(hilbertSeries(extmodule1, Reduce => true), T^(opts.Degree)); -- Here extmodule1 = Ext^{n-1-cohodeg)(M,w).
	-- This is a bounded below series in T. It remains to replace T by U = T^(-1),
	-- compute some properties, and give the output in the form we want. (If HH^cohodeg(X, F(*)) is bounded above and below,
	-- we will write the output in terms of T.)
	if value numerator hilb1 == 0 then (
	    bottomdeg = infinity;
	    topdeg = -infinity;
	    positiveseries = 0;
	    negativeseries = 0)
	else (  -- Here the output is not zero.
	    topdeg = -first min exponents value numerator hilb1; -- This is minus the bottom degree of the numerator,
	    -- which is a Laurent polynomial in T.
	    if value denominator hilb1 == 1 then (  -- In this case, HH^(cohodeg)(X, F(*)) is bounded above and below.
		bottomdeg = -first degree value numerator hilb1;
		positiveseries = substitute(hilb1, T => T^(-1));
		negativeseries = 0)
	    else (
		bottomdeg = -infinity;
		positiveseries = 0;
		negativeseries = substitute(hilb1, T=>U)
		)
	    )
	);
    ("The following is correct in all degrees. Bottom degree:", bottomdeg, "top degree:", topdeg, "cohomology as a series in T:",
	positiveseries, "plus cohomology as a series in U = T^(-1):", negativeseries))
