newPackage(
    "Msolve",
    Version => "1.26.08",
    Date => "August 2026",
    Authors => {
        { Name => "Martin Helmer",  Email => "martin.helmer@swansea.ac.uk", HomePage => "http://martin-helmer.com/" },
        { Name => "Mike Stillman",  Email => "mike@math.cornell.edu",       HomePage => "https://math.cornell.edu/michael-e-stillman" },
        { Name => "Anton Leykin",   Email => "leykin@math.gatech.edu",      HomePage => "https://antonleykin.math.gatech.edu/" },
        { Name => "Mahrud Sayrafi", Email => "mahrud@mcmaster.ca",          HomePage => "https://mahrud.github.io/" },
        },
    Keywords => {"Groebner Basis Algorithms" , "Interfaces"},
    Headline => "interface to the msolve library for solving multivariate polynomial systems using Groebner Bases",
    PackageImports => { "Complexes", "Elimination", "Saturation" },
    AuxiliaryFiles => true,
    DebuggingMode => false
    )

---------------------------------------------------------------------------

export{
    "msolveGB",
    "msolveSyzygy",
    "msolveResolution",
    "MsolveResolution",
    "MsolveDifferential",
    "msolveSaturate",
    "msolveEliminate",
    "msolveRUR",
    "msolveLeadMonomials",
    "msolveRealSolutions",
    "msolveSetup",
    "QQi",
    }

importFrom_Core {
    "raw",
    "rawMatrixReadMsolveFile",
    "rawMatrixWriteMsolveFile",
    --
    "rawMsolveGB",
    "rawMsolveGBRestrictToQuotient",
    "rawMsolveSyzygy",
    "rawMsolveResolution",
    "rawMsolveSaturate",
    "rawMsolveMinimalBetti",
    "rawMsolvePresent",
    "rawResolutionGetFree",
    "rawResolutionGetMatrix",
    "fullgens",
    "degreeToHeft",
    --
    "ContainmentHooks",
    }

-- Direct in-memory interface to msolve's F4, with no temporary files and no
-- string conversion in either direction. Everything the engine cannot take
-- still goes through the executable; these routines return null in that case.

-- mirrors grevlexWeights in e/interface/msolve.cpp: a single GRevLex block over
-- all the variables, with positive weights. A weighted order is fine because
-- the engine applies the substitution x_i -> x_i^(w_i) to the exponents as it
-- marshals them, which carries the order over to msolve's exactly.
isGRevLexEngineRing = R -> (
    mo := toList (options monoid R).MonomialOrder;
    blocks := select(mo, x -> x#0 === GRevLex);
    #blocks === 1
    and #(blocks#0#1) === numgens R
    and all(blocks#0#1, w -> w >= 1)
    and all(mo, x -> member(x#0, {GRevLex, MonomialSize, Position})))

-- The one order msolve implements beyond plain grevlex is a two block degree
-- reverse lexicographic order, and it takes it as the length of the leading
-- block rather than as an order: this is rawMsolveGB's elimination block
-- length, and it makes msolve compare by the degree of the first block, then
-- reverse lexicographically within it, then likewise for the second (the _be
-- comparators in order.c in msolve's neogb).  That is exactly Macaulay2's
-- MonomialOrder => {GRevLex => a, GRevLex => b}: the reduced bases agree
-- element for element, not merely up to a change of order.  So such a ring is
-- usable after all -- not by handing msolve the order, which it cannot
-- represent, but by computing in the plain grevlex ring on the same variables
-- and telling it where the blocks are split.  See msolveBlockGB.  Weights are
-- fine for the same reason as above, the substitution being order preserving
-- one block at a time, and the block degrees msolve accumulates from the
-- scaled exponents being the weighted ones.
--
-- Returns the length of the leading block, or null if the order is not two
-- positively weighted grevlex blocks covering the variables.  Note an
-- Eliminate n order is *not* one of these: Macaulay2 spells it with a weight
-- vector, and the basis it asks for is a different one.
grevlexBlockLength = R -> (
    mo := toList (options monoid R).MonomialOrder;
    blocks := select(mo, x -> x#0 === GRevLex);
    if #blocks === 2
    and #(blocks#0#1) + #(blocks#1#1) === numgens R
    and all(join(blocks#0#1, blocks#1#1), w -> w >= 1)
    and all(mo, x -> member(x#0, {GRevLex, MonomialSize, Position}))
    then #(blocks#0#1))

-- the coefficients msolve can represent: a prime field of positive
-- characteristic below 2^31.  Note isField excludes the coefficients of an
-- unflattened tower such as (ZZ/p[x,y])[u,v]; those still go through
-- toMsolveRing, which flattens them first
msolveCoefficientsUsable = kk -> (
    isField kk
    and char kk != 0
    and char kk <= 2^31
    and precision kk === infinity
    and not instance(kk, GaloisField))

msolveEngineUsable = m -> (
    if not rawMsolvePresent() then return false;
    R := ambient first flattenRing ring m;
    instance(R, PolynomialRing)
    and msolveCoefficientsUsable coefficientRing R
    and isGRevLexEngineRing R)

---------------------------------------------------------------------------

-- A Betti tally key is (level, degree, heft degree), and the degree is a
-- *vector*: the layout above spells the length one case as a single number,
-- and it is only because one number can be subtracted from another that such a
-- table can be laid out slanted and dense, one row per degree minus level.
-- rawMsolveMinimalBetti reports the multidegrees themselves, which have nothing
-- to slant by -- a multidegree minus a level is not a multidegree -- and which
-- are far too sparse to lay out densely, a grading of length seven having its
-- nonzero entries scattered over a lattice.  So what comes back is the nonzero
-- entries and nothing else, grouped by homological level:
--   [r, len, n_0, ..., n_(len-1), (d_1, ..., d_r, heft, beta) x n_i per level i]
-- where r is the degree length, len the homological length, and n_i the number
-- of entries at level i.  This is a plain BettiTally: a table whose degrees are
-- longer than one prints folded onto the heft degree, as BettiTally's own net
-- does, until the caller asks for `multigraded`.
unpackMsolveBetti = w -> (
    (r, len) := (w#0, w#1);
    (stride, offset) := (r + 2, len + 2);
    new BettiTally from flatten for i to len - 1 list (
	recs := pack(stride, take(w, {offset, offset + w#(2+i) * stride - 1}));
	offset = offset + w#(2+i) * stride;
	for rec in recs list (i, take(rec, r), rec#r) => rec#(r+1)))

---------------------------------------------------------------------------

-- A DegreeLimit is a degree of the ring and so may be a multidegree, while
-- msolve schedules its rounds by the single integer heft . deg. Collapsing
-- the one to the other is exactly what Macaulay2 does to its own DegreeLimit
-- before the engine ever sees it, so degreeToHeft is borrowed rather than
-- reimplemented. {} is gb's spelling of "no limit" and 0 is msolve's.
msolveHeftLimit = (R, d) -> (
    h := degreeToHeft(first flattenRing R, d);
    if h === {} then 0 else if first h > 0 then first h
    else error "msolve: expected a degree limit of positive heft")

-- A stopping condition is spelled as msolve's 0-means-no-limit int32.
msolveEngineNat = n -> if n === infinity or n === null then 0 else (
    if instance(n, ZZ) and n >= 0 then n
    else error "msolve: expected a nonnegative integer or infinity")

-- A Groebner basis of the columns of m, a matrix of any number of rows,
-- eliminating the first elim variables. Over R = S/J the generators of J go to
-- msolve as ordinary input placed at component 0, the component its module
-- hash table reserves for the ring: one stored copy of J acts in every
-- component, so nothing is lifted here and nothing is tensored with
-- id_(target m0) once per row. What comes back is already a basis of <m> in
-- R^r, the elements whose lead term lies in in(J) having been dropped inside
-- the engine by the same redundancy marking that drops every other divisible
-- lead term. Compare msolveSaturateEngine, which still lifts.
msolveGBMatrix = (m, elim, opts) -> (
    if not msolveEngineUsable m then return null;
    -- msolve's module F4 offers position over term and term over position on
    -- the components, but not an induced (Schreyer) order: res.h's
    -- RES_MORD_SCHREYER wants per component base monomials that only its
    -- resolution engine supplies, and export_module_f4 has no way to accept
    -- them.  What it would return under a Schreyer target still generates the
    -- right submodule, but with the lead terms of the wrong order, and
    -- msolveDeclareGB's forceGB would declare a non-Groebner basis to be one,
    -- silently corrupting poincare, %, and everything else lead term based.
    -- Twists of the target are fine: rawMsolveModuleGB passes no row degrees,
    -- so its order compares deg(m) and then the component, as Macaulay2 does.
    if numrows m > 1 and schreyerOrder target m != 0 then return null;
    threads := opts.Threads ?? allowableThreads;
    verbosity := opts.Verbosity ?? gbTrace;
    -- callers that have no DegreeLimit in their option table ask for the whole
    -- basis, which is what every one of them did before the option existed
    deglimit := if opts.?DegreeLimit then msolveHeftLimit(ring m, opts.DegreeLimit) else 0;
    (R0, phi) := flattenRing ring m;
    G := rawMsolveGB(raw sub(matrix m, R0), elim, deglimit, threads, verbosity);
    -- a flat ring needs neither of the two maps back, and asking for them
    -- would cost a copy of the basis apiece
    if R0 === ring m then return map(target m, , G);
    map(target m, , phi map(R0, G)))

msolveGBEngine = (m, elim, opts) -> msolveGBMatrix(m, elim, opts)

-- F4SAT additionally needs the characteristic to be larger than 2^16, see
-- https://github.com/algebraic-solving/msolve/issues/165
--
-- Unlike msolveGBMatrix this still lifts a quotient rather than handing J to
-- the engine: F4SAT is the ideal engine, whose hash table has no component
-- slot and so no component 0 for a ring element to live in, and it wants
-- exactly the polynomials to saturate by as its trailing generators. Lifting
-- is correct because the preimage in S of saturate(I*R0, f) is (I + J) : f^oo,
-- so J goes in with the generators and the columns whose lead term already
-- lies in in(J) come back out again afterwards.
msolveSaturateEngine = (m, f, opts) -> (
    if not msolveEngineUsable m or char ring m < 2^16 then return null;
    (R0, phi) := flattenRing ring m;
    S := ambient R0;
    threads := opts.Threads ?? allowableThreads;
    verbosity := opts.Verbosity ?? gbTrace;
    if S === R0 then return map(S, rawMsolveSaturate(
            raw sub(matrix m, S), raw sub(matrix {{f}}, S), threads, verbosity));
    m0 := lift(sub(matrix m, R0), S);
    f0 := lift(sub(matrix {{f}}, R0), S);
    G := rawMsolveSaturate(raw(m0 | presentation R0), raw f0, threads, verbosity);
    -- keep only the columns not already accounted for by in(J); a column is
    -- accounted for exactly when its lead term is already in R0's own initial
    -- ideal.  R0 already has that as a MonomialTable, built when it was
    -- constructed for its own ring arithmetic, so this is a membership query
    -- against existing engine data rather than a second Groebner basis
    -- computation of leadTerm(presentation R0) done at the interpreter level.
    map(target m, , phi substitute(
            map(S, rawMsolveGBRestrictToQuotient(G, raw R0)), vars R0)))

-- used in msolveEliminate
importFrom_Core { "monoidIndices" }
importFrom_Elimination { "eliminationRing" }

---------------------------------------------------------------------------

msolveMinimumVersion = "0.7.0"
msolveProgram = findProgram("msolve", "msolve -h",
    MinimumVersion => (msolveMinimumVersion, "msolve -V"),
    RaiseError => false,
    Verbose => debugLevel > 0)

if msolveProgram === null then (
    printerr("warning: could not find msolve with version at least v" | msolveMinimumVersion);
    -- note: msolve -h returns status code 1 :/
    msolveProgram = findProgram("msolve", "true", Verbose => debugLevel > 0))

msolveDefaultOptions = new OptionTable from {
    Threads => null, Verbosity => null }

---------------------------------------------------------------------------
-- Which of Macaulay2's gb, syz and res options msolve can honor
---------------------------------------------------------------------------
-- msolveDefaultOptions above is deliberately tiny: Threads and Verbosity are
-- the only two knobs every msolve entry point, executable or engine, already
-- takes. The table below records what it would take to accept the rest, so
-- that the next person to reach for one of these knows whether it is a matter
-- of forwarding an argument or of changing msolve itself.
--
-- The engine entry points are export_module_f4, export_module_frame,
-- export_module_resolution, export_module_betti and the res_comp_t handle in
-- msolve's src/neogb/res.h, reached from here through rawMsolveGB,
-- rawMsolveSyzygy, rawMsolveMinimalBetti, rawMsolvePoincare and
-- rawMsolveResolution. Between them they accept max_level, syz_of, minimal,
-- verify, ht_size, nr_threads, max_nr_pairs, la_option, reduce_gb and
-- info_level, plus the res_stop_t of the entry points that take one -- and
-- nothing else. res_stop_t holds the three stops msolve has: a degree
-- ceiling, a syzygy count and a row bound. Every other "stop after N of
-- something" below is still engine work rather than plumbing.
--
-- None of the three makes a computation resumable, and that is not an
-- oversight to be plumbed around: finalize_f4 frees the pair set and the
-- hash table, and reduce_final_basis rebuilds and renumbers the basis, so
-- the reduced basis msolve returns and the basis a continuation would need
-- are different objects. A truncated answer is therefore never cached as a
-- complete one -- see msolveGBHook, which still refuses DegreeLimit for
-- exactly that reason.
--
--   Done      forwarded already; see the option it is spelled as
--   Yes       already reachable; forwarding it is interface work here
--   Partial   an msolve knob exists but does not mean the same thing
--   Engine    wants a change inside msolve; see the milestone in the plan
--   No        M2-side or meaningless for F4; nothing to forward
--
-- gb and syz:
--
--   Option                     |         | Notes
--   ---------------------------|---------|--------------------------------------
--   Algorithm                  | No      | msolve is F4 only; la_option picks the
--                              |         | linear algebra variant, not the pair
--                              |         | selection strategy
--   BasisElementLimit          | Engine  | no stop conditions
--   ChangeMatrix               | Engine  | the graph module (f_j, e_j) already
--                              |         | carries it: the adjoined part of each
--                              |         | GB element is its expression in the
--                              |         | input. RES_SYZ_OF_INPUT computes this
--                              |         | GB and discards exactly the elements a
--                              |         | change matrix would keep
--   CodimensionLimit           | Engine  | no stop conditions
--   DegreeLimit                | Done    | res_stop_t's max_degree, a ceiling on
--                              |         | the F4 round loop, which was already
--                              |         | degree by degree (symbol.c selects the
--                              |         | minimal degree). It is a multidegree,
--                              |         | coarsened to its heft because that is
--                              |         | what the schedule counts up; spelled as
--                              |         | msolveGB's DegreeLimit. Not resumable,
--                              |         | so not offered to the gb hook
--   GBDegrees                  | Partial | isGRevLexEngineRing passes a weight
--                              |         | vector, applied as the substitution
--                              |         | x_i -> x_i^(w_i). The Betti and
--                              |         | Hilbert entry points no longer do:
--                              |         | they hand msolve the ring's grading
--                              |         | outright -- a degree for each
--                              |         | variable, a heft, and any torsion
--                              |         | factors -- and it computes in the
--                              |         | heft order that induces, which is
--                              |         | what lets a multigraded ring through.
--                              |         | The Groebner basis paths still
--                              |         | substitute; see collectGrading in
--                              |         | e/interface/msolve.cpp
--   HardDegreeLimit            | Engine  | as DegreeLimit, plus discarding above it
--   Hilbert                    | Engine  | msolve has no Hilbert driven early
--                              |         | termination. Now that M5 *produces* the
--                              |         | numerator (rawMsolvePoincare), teaching
--                              |         | the F4 to *consume* one is the natural
--                              |         | next step and the highest value item here
--   MaxReductionCount          | No      | an M2 gb auto-reduction knob; F4 reduces
--                              |         | a whole Macaulay matrix at once
--   PairLimit                  | Partial | max_nr_pairs exists but is st->mnsel,
--                              |         | the most pairs selected in one round
--                              |         | (symbol.c:268), a batching knob. It is
--                              |         | not a stop and must not be mapped to one
--   StopBeforeComputation      | No      | interpreter side
--   StopWithMinimalGenerators  | Engine  | no stop conditions
--   Strategy                   | Partial | la_option and reduce_gb are the only
--                              |         | strategy dials msolve exposes
--   SubringLimit               | Engine  | msolve takes an elimination block length
--                              |         | but has no count to stop at
--   Syzygy                     | Yes     | export_module_resolution with
--                              |         | syz_of = RES_SYZ_OF_INPUT and
--                              |         | max_level = 2; this is msolveSyzygy,
--                              |         | just not surfaced as a gb option
--   SyzygyLimit                | Done    | res_stop_t's syz_limit, spelled as
--                              |         | msolveSyzygy's SyzygyLimit. A genuine
--                              |         | early stop rather than a cap on the
--                              |         | output: under position over term the
--                              |         | lead term alone says whether a basis
--                              |         | element is a syzygy, so the round loop
--                              |         | counts them for one comparison apiece.
--                              |         | A round is always finished, so the
--                              |         | extras of the last one are discarded
--   SyzygyRows                 | Done    | res_stop_t's syz_rows, spelled as
--                              |         | msolveSyzygy's SyzygyRows. Free as an
--                              |         | output filter, since the rows are the
--                              |         | adjoined components; pruning the *work*
--                              |         | still needs those components ordered
--                              |         | last, which is the elimination block
--                              |         | question the module orders do not yet
--                              |         | answer. The result is a submatrix of the
--                              |         | syzygy matrix, so M times it is not zero
--
-- res:
--
--   DegreeLimit                | Yes     | res_stop_t's max_degree already reaches
--                              |         | the Groebner basis the frame is built
--                              |         | on, and res_diff_compute is already
--                              |         | driven one (level, degree) at a time;
--                              |         | what is missing is a res_stop_t on
--                              |         | export_module_betti and res_comp_new,
--                              |         | plus deciding what a truncated Betti
--                              |         | table means for the invariants
--   HardDegreeLimit            | Engine  | as above
--   LengthLimit                | Done    | max_level, taken by every entry point
--                              |         | that builds a frame; plumbed as
--                              |         | rawMsolveMinimalBetti's length_limit and
--                              |         | as msolveResolution's LengthLimit.
--                              |         | Note whole module invariants (poincare,
--                              |         | pdim, regularity, dim, degree) are refused
--                              |         | on a truncated table, since the frame is
--                              |         | nonminimal and can run past numgens R
--   PairLimit                  | Partial | as above
--   ParallelizeByDegree        | No      | msolve parallelizes inside each linear
--                              |         | algebra step instead; use Threads
--   SortStrategy               | No      | the block order within a frame level is
--                              |         | fixed at degree ascending then ring order
--                              |         | descending, which is what M2 does anyway
--   StopBeforeComputation      | No      | interpreter side
--   Strategy                   | Yes     | export_module_betti's minimal flag is
--                              |         | Minimal vs Nonminimal: with minimal = 0
--                              |         | it reports frame ranks and performs no
--                              |         | field arithmetic past the Groebner basis
--   SyzygyLimit                | Partial | res_stop_t's syz_limit counts syzygies
--                              |         | of the *input*, which is level 2 only;
--                              |         | res's SyzygyLimit is a bound at every
--                              |         | level and the frame has no counter yet
--
-- msolve also offers three knobs M2's option tables have no name for:
-- ht_size (initial hash table size, a log2), verify (an exact d o d = 0 check
-- over the whole complex, several times the cost of the resolution itself),
-- and reduce_gb (whether to fully interreduce the basis before returning).

runMsolve = (mIn, mOut, args, opts) -> runProgram(msolveProgram,
    demark_" " { args,
	"-t", toString(opts.Threads ?? allowableThreads),
	"-v", toString(opts.Verbosity ?? gbTrace),
	"-f", toString mIn,
	"-o", toString mOut },
    KeepFiles => true,
    RaiseError => true,
    Verbose => (opts.Verbosity ?? gbTrace) > 0)

-- e.g. turns x_(0,0)... to p_0...
toMsolveRing = I -> (
    S0 := ring I;
    S1 := first flattenRing S0;
    kk := ultimate(coefficientRing, S0);
    if not instance(S1, PolynomialRing) or instance(kk, GaloisField)
    or not isField kk or char kk > 2^31 or precision kk < infinity
    then error "msolve: expected an ideal in a polynomial ring over QQ or ZZ/p with characteristic less than 2^31";
    -- resets the variables to p_0...
    S := newRing(S1, Variables => numgens S1);
    S, kk, substitute(I, vars S))

toMsolveString = X -> (
    elts := if instance(X, List)     then X
    else if instance(X, Ring)        then X_*
    else if instance(X, Ideal)       then X_*
    else if instance(X, RingElement) then {X};
    str := toExternalString elts;
    -- (2/5)*x -> 2/5*x
    str = replace("[)(]", "", str);
    -- {x,y,z} -> x,y,z
    str_(1, #str-2))

toMsolveInput = (S, K, I) -> demark_newline {
    toMsolveString S,
    toString char K,
    replace(",", ",\n",
	toMsolveString I)}

-- over ZZ/p the engine can write the input file directly from the matrix,
-- which avoids building a (potentially enormous) string in the interpreter;
-- over QQ the engine has no rational coefficients, so we fall back to strings.
use'writeMsolveInputFile := true;
writeMsolveInputFile = (S, K, I, mIn) -> (
    if use'writeMsolveInputFile and char K > 0 and #I > 0
    then rawMatrixWriteMsolveFile(raw matrix {I}, mIn)
    else (mIn << toMsolveInput(S, K, I) << endl << close;);
    mIn)

msolve = (S, K, I, args, opts) -> (
    tmp := temporaryFileName();
    mIn := writeMsolveInputFile(S, K, I, tmp | "-in.ms");
    mOut := tmp | "-out.ms";
    runMsolve(mIn, mOut, args, opts);
    mOut)

use'readMsolveOutputFile := true;
readMsolveOutputFile = method()
readMsolveOutputFile(Ring,String) := Matrix => (R,mOut) -> if use'readMsolveOutputFile 
    -- TODO: this substitution should be unnecessary,
    -- but without it the result for tower rings is in the wrong order!
    then substitute(map(R, rawMatrixReadMsolveFile(raw R, mOut)), R) else (
	--the line below should be replaced by a call to the C-function to parse the string
	-- this is a hack that has global consequences (e.g. breaks rings with p_i vars)
	use newRing(R, Variables => numgens R);
	substitute(matrix {value readMsolveList get mOut}, vars R))

readMsolveList = mOutStr -> (
    mOutStr = toString stack select(lines mOutStr,
	line -> not match("#", line));
    mOutStr = replace("\\[", "{", mOutStr);
    mOutStr = replace("\\]", "}", mOutStr);
    -- e.g. 'p_0' to "p_0"
    mOutStr = replace("'", "\"",  mOutStr);
    mOutStr = first separate(":", mOutStr);
    mOutStr)

---------------------------------------------------------------------------
-- Core msolve algorithm calls
---------------------------------------------------------------------------

-- DegreeLimit is not resumable: msolve keeps no state between calls, and its
-- reduced basis is a different object from the one a continuation would need
-- (reduce_final_basis renumbers the whole basis and free_meta_data drops the
-- pair set). Asking again for a larger limit recomputes from the input, and a
-- truncated basis is never cached as though it were complete -- which is why
-- msolveGBHook below still refuses the option, and why a limited msolveGB
-- Ideal returns the matrix rather than forcing a GroebnerBasis object.
msolveGBOptions = new OptionTable from {
    Threads => null, Verbosity => null, DegreeLimit => {} }

msolveGB = method(TypicalValue => Matrix, Options => msolveGBOptions)
msolveGB Matrix := opts -> M -> (
    if (G := msolveGBEngine(M, 0, opts)) =!= null then return G;
    error "msolveGB: expected a matrix over a GRevLex polynomial ring over ZZ/p, with 0 < p < 2^31")
msolveGB Module := opts -> M -> (
    if M.?relations
    then msolveGB(fullgens M, opts) -- TODO: do we need SyzygyRows => numgens source generators M?
    else msolveGB(generators M, opts))
msolveGB Ideal := opts -> I0 -> (
    if (G := msolveGBEngine(generators I0, 0, opts)) =!= null
    then return if opts.DegreeLimit === {} then gens forceGB G else G;
    --
    if opts.DegreeLimit =!= {} then error(
        "msolveGB: a degree limit needs the msolve library; ",
        "the msolve executable computes a complete basis and nothing less");
    (S, K, I) := toMsolveRing I0;
    mOut := msolve(S, K, I_*, "-g 2", opts);
    gens forceGB readMsolveOutputFile(ring I0, mOut))

-- The syzygies come from msolve's module F4, which takes a quotient ring but
-- not a block order -- compare msolveBlockGB, which takes the block order but
-- only for one row.  Unlike msolveEngineUsable this does not flatten, nothing
-- below mapping the result back, so a tower is left to Macaulay2; the ambient
-- of a tower is a polynomial ring over a polynomial ring and so is refused by
-- msolveCoefficientsUsable anyway, as is ZZ/p, whose ambient is ZZ.
msolveSyzygyUsable = m -> (
    rawMsolvePresent()
    and (instance(R := ring m, PolynomialRing) or instance(R, QuotientRing))
    and instance(A := ambient R, PolynomialRing)
    and msolveCoefficientsUsable coefficientRing A
    and isGRevLexEngineRing A)

-- SyzygyRows keeps only the first rows of the syzygy matrix, so the result is
-- a submatrix and M * msolveSyzygy(M, SyzygyRows => n) is not zero; the
-- dropped rows come back as zero rows, so the target is unchanged. This is
-- what a caller wants when the columns of M are generators followed by
-- relations and only the coefficients on the generators matter, which is the
-- SyzygyRows => numgens source generators M of gb.m2.
msolveSyzygyOptions = new OptionTable from {
    Threads => null, Verbosity => null,
    DegreeLimit => {}, SyzygyLimit => infinity, SyzygyRows => infinity }

msolveSyzygy = method(TypicalValue => Matrix, Options => msolveSyzygyOptions)
msolveSyzygy Matrix := opts -> M -> (
    if not rawMsolvePresent() then
        error "msolveSyzygy: this Macaulay2 was built without the msolve library";
    if not msolveSyzygyUsable M
    then error "msolveSyzygy: expected a matrix over a GRevLex polynomial ring over ZZ/p, with 0 < p < 2^31";
    R := ring M;
    threads := opts.Threads ?? allowableThreads;
    verbosity := opts.Verbosity ?? gbTrace;
    -- msolve spells "all of them" as 0 and Macaulay2 spells it as infinity,
    -- so an explicit 0 -- gb.m2's gbOnly asks for exactly that -- has to be
    -- answered here: no rows kept means every column is dropped
    if opts.SyzygyRows === 0 then return map(source M, R^0, 0);
    -- a bound at or past the number of columns of M is no bound at all
    syzrows := msolveEngineNat opts.SyzygyRows;
    if syzrows >= numcols M then syzrows = 0;
    map(R, rawMsolveSyzygy(raw matrix M,
            msolveEngineNat opts.SyzygyLimit, syzrows,
            msolveHeftLimit(R, opts.DegreeLimit), threads, verbosity)))

---------------------------------------------------------------------------
-- Free resolutions, one free module and one differential at a time
---------------------------------------------------------------------------

-- msolveResolution returns a live computation rather than a complex.  Building
-- it runs the module Groebner basis and the whole Schreyer frame; the frame is
-- combinatorial, so from that point on C_i is free -- rank and degrees, no
-- field arithmetic -- while C.dd_i is what makes msolve reduce, and then only
-- up to level i.  Levels already computed are remembered, so asking twice
-- costs nothing and asking out of order costs no more than asking in order.
--
-- The complex is the *nonminimal* one, as with res(..., Strategy => Nonminimal):
-- C_1 is the Groebner basis of the image of the input, not the input columns,
-- since msolve keeps no change of basis between them.  Its ranks therefore
-- depend on the Groebner basis and so on the module order, and need not agree
-- with Macaulay2's own nonminimal resolution of the same module -- only the
-- minimal Betti numbers are an invariant, and those come from
-- rawMsolveMinimalBetti without any of this being materialized.

protect RawComputation

MsolveResolution = new Type of MutableHashTable
MsolveResolution.synonym = "msolve resolution"

MsolveDifferential = new Type of MutableHashTable
MsolveDifferential.synonym = "differential of an msolve resolution"

msolveResolutionOptions = new OptionTable from {
    Threads => null, Verbosity => null, LengthLimit => infinity }

msolveResolution = method(TypicalValue => MsolveResolution,
    Options => msolveResolutionOptions)
msolveResolution Matrix := opts -> M -> (
    if not rawMsolvePresent() then
        error "msolveResolution: this Macaulay2 was built without the msolve library";
    R := ring M;
    if not instance(R, PolynomialRing)
    or not isField(coefficientRing R)
    or char R == 0 or char R >= 2^31
    or precision(coefficientRing R) =!= infinity
    or instance(coefficientRing R, GaloisField)
    or not isGRevLexEngineRing R
    then error "msolveResolution: expected a matrix over a GRevLex polynomial ring over ZZ/p, with 0 < p < 2^31";
    len := if opts.LengthLimit === infinity then 0 else opts.LengthLimit;
    if not instance(len, ZZ) or len < 0
    then error "msolveResolution: expected LengthLimit to be a nonnegative integer or infinity";
    threads := opts.Threads ?? allowableThreads;
    verbosity := opts.Verbosity ?? gbTrace;
    G := rawMsolveResolution(raw matrix M, len, threads, verbosity);
    -- msolve reports the reason on stderr; the overwhelmingly likely one is
    -- that a resolution is a graded object and this input is not homogeneous
    if G === null then error("msolveResolution: msolve could not resolve this "
        | "matrix; it must be homogeneous over a singly graded ring whose "
        | "variable degrees are its GRevLex weights");
    C := new MsolveResolution from {
        symbol ring => R,
        symbol source => M,
        RawComputation => G,
        symbol cache => new CacheTable };
    C.dd = new MsolveDifferential from { symbol target => C };
    C)
msolveResolution Ideal  := opts -> I -> msolveResolution(generators I, opts)
-- the resolution is of the cokernel, so it is the presentation that goes in
msolveResolution Module := opts -> N -> msolveResolution(presentation N, opts)

ring MsolveResolution := C -> C.ring

-- F_i, free and immediate: this is the frame
MsolveResolution _ ZZ := Module => (C, i) -> (
    if i < 0 then return (ring C)^0;
    if C.cache#?i then return C.cache#i;
    C.cache#i = new Module from (ring C, rawResolutionGetFree(C.RawComputation, i)))

-- d_i, and the only thing here that costs anything
MsolveDifferential _ ZZ := Matrix => (D, i) -> (
    C := D.target;
    if i < 1 or C_i == 0 then return map(C_(i-1), C_i, 0);
    map(C_(i-1), C_i, rawResolutionGetMatrix(C.RawComputation, i)))

-- the last level with a nonzero free module; the frame is nonminimal, so this
-- is not bounded by numgens R, and with a LengthLimit it is the limit itself
length MsolveResolution := C -> (
    if C.cache#?(symbol length) then return C.cache#(symbol length);
    i := 0;
    while C_(i+1) != 0 do i = i+1;
    C.cache#(symbol length) = i)

net MsolveResolution := C -> horizontalJoin between(" <-- ",
    apply(1 + length C, i -> net C_i))

-- The whole thing at once, as an ordinary Complex.  This materializes every
-- differential, which is precisely what the live handle exists to avoid, so it
-- is something the caller asks for rather than something another method does
-- behind its back -- C_i and length C stay free.
--
-- What it buys is everything the Complexes package knows how to do with a free
-- resolution, and in particular `minimize`, which turns the nonminimal
-- resolution msolve computes into the minimal one.  The result resolves
-- coker of the input matrix, so `betti minimize complex C` is the minimal
-- Betti table -- the same one rawMsolveMinimalBetti gets from ranks alone,
-- without materializing anything, and far more cheaply.
complex MsolveResolution := Complex => {Base => 0} >> opts -> C -> (
    n := length C;
    if n === 0 then complex(C_0, Base => opts.Base)
    else complex(apply(toList(1 .. n), i -> C.dd_i), Base => opts.Base))

importFrom_Core "numallvars"
msolveLeadMonomials = method(TypicalValue => Matrix, Options => msolveDefaultOptions)
msolveLeadMonomials Ideal := opts -> I0 -> (
    -- TODO: premute the coefficient variables to make this work for tower rings
    S0 := ring I0;
    if numgens S0 =!= S0.numallvars then error "msolveLeadMonomials: unsupported tower ring";
    -- msolve's -g 1 computes the same basis and merely prints less of it, so
    -- with the engine we take the lead terms of the basis it hands back
    if (G := msolveGBEngine(generators I0, 0, opts)) =!= null
    then return gens forceGB leadTerm G;
    (S, K, I) := toMsolveRing I0;
    mOut := msolve(S, K, I_*, "-g 1", opts);
    gens forceGB readMsolveOutputFile(ring I0, mOut))
-- TODO: add leadMonomials Ideal, then add this as a hook

msolveEliminate = method(Options => msolveDefaultOptions)
msolveEliminate(RingElement, Ideal) := Ideal => opts -> (elimvar,  I) -> msolveEliminate(I, {elimvar}, opts)
msolveEliminate(List,        Ideal) := Ideal => opts -> (elimvars, I) -> msolveEliminate(I,  elimvars, opts)
msolveEliminate(Ideal, RingElement) := Ideal => opts -> (I,  elimvar) -> msolveEliminate(I, {elimvar}, opts)
msolveEliminate(Ideal,        List) := Ideal => opts -> (I0, elimvars) -> (
    -- turns generators into indices, but also allows indices
    elimIndices := monoidIndices(S0 := ring I0, elimvars);
    keepvars := S0_* - set S0_*_elimIndices;
    -- gives ring maps to and from a ring with elimvars first, keepvars last
    (toS0', toS0) := eliminationRing(elimIndices, S0);
    (S, K, I) := toMsolveRing(I' := toS0' I0);
    S' := if char K === 0
    then K(monoid [keepvars]) -- msolve does not return the remaining generators over QQ
    else newRing(ring I', MonomialOrder => {#elimvars, #keepvars}); -- but over ZZ/p it does
    -- toMsolveRing already moved the ideal into a plain grevlex ring, which is
    -- what msolve wants: the elimination is expressed by the block length, not
    -- by the monomial order of the ring
    if (G := msolveGBEngine(generators I, #elimIndices, opts)) =!= null
    then return ideal substitute(G, vars S');
    mOut := msolve(S, K, I_*, "-g 2 -e " | length elimIndices, opts);
    ideal readMsolveOutputFile(S', mOut))

msolveSaturate = method(TypicalValue => Matrix, Options => msolveDefaultOptions)
msolveSaturate(Ideal, RingElement) := opts -> (I0, f0) -> (
    if (G := msolveSaturateEngine(generators I0, f0, opts)) =!= null
    then return gens forceGB G;
    (S, K, I) := toMsolveRing I0;
    f := substitute(f0, vars S);
    -- see https://github.com/algebraic-solving/msolve/issues/165
    if char K < 2^16 or 2^31 < char K then error "msolveSaturate: expected characteristic between 2^16 and 2^31 for F4SAT";
    -- msolve expects a list of the generators of the ideal followed by f
    mOut := msolve(S, K, I_* | {f}, "-S -g 2", opts);
    gens forceGB readMsolveOutputFile(ring I0, mOut))

--------------------------------------------------------------------------------
-- Routing Macaulay2's own commands through msolve
--------------------------------------------------------------------------------

-- Everything below is opt in: calling msolveSetup() redirects gb,
-- groebnerBasis, eliminate and kernel to msolve whenever msolve can handle the
-- input, and leaves Macaulay2's own implementation in charge otherwise.
-- Compare with fast-kernel.m2, which this generalizes.

M2DefaultGB        = lookup(gb, Matrix)
M2DefaultGBasis    = lookup(groebnerBasis, Matrix)
M2DefaultEliminate = lookup(eliminate, List, Ideal)
M2DefaultMinimalBetti = lookup(minimalBetti, Module)

-- true if msolve is applicable to the one-rowed matrix m
msolveApplicable = m -> m != 0 and msolveEngineUsable m

-- Declares the columns of G a Groebner basis.  Note forceGB caches the
-- declaration on G itself, whereas `generators forceGB G` rebuilds the matrix
-- from the raw computation and hands back one with an empty cache: it is G
-- that has to be kept if a later gb of it is to cost nothing.
msolveDeclareGB = G -> (
    G = forceGB G;
    -- records that no Hilbert function hint was needed, as gb.m2 would
    G#"rawGBSetHilbertFunction log" = true;
    G)

-- A GroebnerBasis object for m, computed by msolve and declared with forceGB.
-- Returns null when msolve does not apply, so callers can fall back.
msolveForceGB = (m, elim, opts) -> (
    if not msolveApplicable m then return null;
    if (G := msolveGBMatrix(m, elim, opts)) === null then return null;
    msolveDeclareGB G)

-- The same variables in the same order, but with a plain grevlex order: msolve
-- expresses an elimination by the length of a leading block of variables rather
-- than by the monomial order of the ring, so the block order that
-- eliminationRing produces is not what should be handed to it.
msolveGRevLexRing = R -> (
    A := ambient R;
    (coefficientRing A)(monoid [Variables => numgens A, Degrees => degrees A]))

-- A Groebner basis of the one-rowed matrix m for the order eliminating the
-- first e variables, returned as a pair (S1, G) with G a matrix over the plain
-- grevlex ring S1 on the same variables. This is the trick fast-kernel.m2 uses:
-- the ring m comes from -- whether it is the block ordered ring produced by
-- eliminationRing or by graphIdeal -- carries an order msolve cannot represent,
-- but msolve does not need it, since it takes the elimination as a block length.
-- A quotient contributes its relations to the generators here rather than
-- going to the engine as msolveGBMatrix does, since res.h refuses to combine
-- an elimination block with the module order that component 0 lives in; the
-- relations map to zero on the way back, so nothing filters them out.
-- Returns null when msolve does not apply.
msolveElimGB = (m, e) -> (
    if not rawMsolvePresent() then return null;
    R1 := ring m;
    A1 := ambient R1;
    S1 := msolveGRevLexRing R1;
    mm := lift(matrix m, A1);
    if A1 =!= R1 then mm = mm | presentation R1;
    m1 := substitute(mm, vars S1);
    if not msolveApplicable m1 then return null;
    G := msolveGBMatrix(m1, e, msolveDefaultOptions);
    if G === null then return null;
    -- an element lies in the elimination ideal exactly when it does not involve
    -- any of the first e variables; those coming from the relations of a
    -- quotient map to zero later on and so contribute nothing
    (S1, matrix(S1, {select(first entries G, f -> all(e, i -> degree(S1_i, f) === 0))})))

-- I intersected with the subring on the variables other than v, computed by
-- msolve and returned as an ideal of ring I. Returns null when msolve does not
-- apply. Compare with eliminate in fast-kernel.m2.
msolveEliminationGB = (I, v) -> (
    R := ring I;
    if I == 0 then return null;
    idx := unique monoidIndices_R v;
    -- puts the variables to be eliminated first, which is the form msolve wants
    (toR1, toR) := eliminationRing(idx, R);
    J := toR1 I;
    result := msolveElimGB(generators J, #idx);
    if result === null then return null;
    (S1, G) := result;
    ideal toR substitute(G, vars ring J))

-- A Groebner basis of the one rowed matrix m for the two block grevlex order of
-- its own ring, computed by msolve and declared with forceGB.  Nothing is
-- approximated: msolve's block order *is* that order, as grevlexBlockLength
-- explains, so what comes back is the reduced basis Macaulay2 would have
-- computed, and a quotient ring built from it reduces correctly.  This is the
-- same detour as msolveElimGB, which computes in msolveGRevLexRing and passes
-- the block length, only kept as a basis of the block ring rather than filtered
-- down to the elimination ideal.
--
-- This is where a tower arrives, and not by accident: a tower is *stored* with
-- a block order, kk[x,y,z][a,b,c] flattening to kk[a,b,c,x,y,z] with the outer
-- variables their own leading grevlex block, and Macaulay2's own Groebner
-- bases over a tower are computed in that very order.  So the block length
-- msolve wants is the number of outer variables, and the only extra work is
-- carrying the matrix across the flattening and its inverse -- which
-- flattenRing hands over, its two maps being mutually inverse isomorphisms.
--
-- Two restrictions remain.  Modules are excluded: msolve does not combine a
-- block order with a module order, and rawMsolveGB refuses the combination.
-- And the ring must be commutative, the plain grevlex ring computed in keeping
-- neither a Weyl nor a skew structure.
--
-- Returns null when msolve does not apply, so callers can fall back.
msolveBlockGB = (m, opts) -> (
    if not rawMsolvePresent() or m == 0 or numrows m =!= 1 then return null;
    if not isCommutative ring m then return null;
    -- phi and psi are the flattening and its inverse; both are the identity
    -- unless the ring is a tower
    (R0, phi, psi) := flattenRing(ring m, Result => (Ring, RingMap, RingMap));
    A := ambient R0;
    -- checked before the plain grevlex ring is built below, so that a gb over a
    -- ring msolve has no use for does not construct one on every call
    if not instance(A, PolynomialRing)
    or not msolveCoefficientsUsable coefficientRing A then return null;
    if (a := grevlexBlockLength A) === null then return null;
    -- the relations of a quotient are appended to the generators, and the whole
    -- thing moved to the plain grevlex ring on the same variables; that ring is
    -- not a quotient, so msolveGBMatrix below hands the engine no component 0
    -- and this path keeps the lift, res.h refusing to combine a block length
    -- with the module order.  msolveEngineUsable is checked there, of that ring
    rels := if A === R0 then null else presentation R0;
    m0 := lift(phi matrix m, A);
    if rels =!= null then m0 = m0 | rels;
    G := msolveGBMatrix(substitute(m0, vars msolveGRevLexRing R0), a, opts);
    if G === null then return null;
    -- back in the block ring, whose order is the one msolve computed in
    G = substitute(G, vars A);
    -- keep only what in(J) does not already account for, exactly as in
    -- msolveSaturateEngine: a membership query against R0's existing
    -- MonomialTable, rather than a second Groebner basis computation of
    -- leadTerm(rels) here
    if rels =!= null then G = map(A, rawMsolveGBRestrictToQuotient(raw G, raw R0));
    msolveDeclareGB map(target m, , psi substitute(G, vars R0)))

msolveGBHook = options gb >> opts -> m -> (
    -- msolve computes a full reduced basis, so a subring or degree limited
    -- request has to go back to Macaulay2's own implementation
    if opts.DegreeLimit =!= {} or opts.SubringLimit =!= infinity
    or opts.ChangeMatrix or opts.Syzygies
    or opts.StopWithMinimalGenerators then return null;
    -- a ring with two grevlex blocks is msolve's own block order, which it takes
    -- as a block length rather than as an order; the two are mutually exclusive
    msolveForceGB(m, 0, msolveDefaultOptions)
    ?? msolveBlockGB(m, msolveDefaultOptions))

-- The syzygies of m, computed by msolve's module F4.
--
-- msolve computes the whole syzygy module in one go and has no stop conditions
-- at all (see the option table above), so the two options that merely shape the
-- answer are honored by trimming what it returns, while the ones that would
-- have to stop the computation early -- a degree, a basis element or a pair
-- count -- send the request back to Macaulay2 rather than being ignored.
--
-- Only homogeneous input is taken: rawMsolveSyzygy reads the syzygies off a
-- graded resolution of one step, which msolve will not start otherwise.  That
-- is also why what comes back is only a generating set and not the basis
-- Macaulay2 would return, in a different order than it would return it -- see
-- the note in rawMsolveSyzygy -- so, as in the Default hook of syz, the
-- generators are minimalized, and then sorted the way syz sorts a basis it
-- reads off a Groebner basis, ascending by degree and then by lead monomial.
-- Minimalizing on its own does not do that: it inherits the order of what it
-- was handed, and what msolve hands over is sorted neither way.
--
-- Returns null when msolve does not apply, so the remaining hooks can be tried;
-- takes its arguments the way runHooks hands them over, as one pair.
msolveSyzygyHook = (opts, m) -> (
    if opts.DegreeLimit =!= {} or opts.HardDegreeLimit =!= null
    or opts.BasisElementLimit =!= infinity or opts.PairLimit =!= infinity
    or opts.StopBeforeComputation
    or m == 0 or not isHomogeneous m or not msolveSyzygyUsable m then return null;
    S := msolveSyzygy m;
    -- SyzygyRows asks for the coefficients on the first few generators only,
    -- which is the image of the syzygy module under the projection onto them;
    -- msolve's columns generate the whole module, so their first rows generate
    -- that, and minimalizing below brings them down to a basis of it
    if opts.SyzygyRows =!= infinity
    then S = S^(toList(0 .. min(opts.SyzygyRows, numrows S) - 1));
    S = if numrows S > 0 then mingens image S else S_(toList(0 .. -1));
    S = sort(S, DegreeOrder => Ascending);
    -- SyzygyLimit asks Macaulay2 to stop once it has that many syzygies, which
    -- leaves it holding the ones of lowest degree.  msolve cannot stop, so that
    -- many of the lowest degree ones stand in for them: a minimal generating
    -- set of the whole module cut down, rather than the partial computation
    -- Macaulay2 would have interrupted
    if opts.SyzygyLimit =!= infinity
    then S = S_(toList(0 .. min(opts.SyzygyLimit, numcols S) - 1));
    S)

-- rawMsolveMinimalBetti only takes a length limit, not a degree limit, and
-- always resolves in the ring's own heft grading (see the option table above),
-- so a DegreeLimit request has to go back to Macaulay2's own implementation;
-- LengthLimit => infinity is spelled 0, meaning the whole resolution
msolveMinimalBettiHook = options minimalBetti >> opts -> M -> (
    if opts.DegreeLimit =!= null then return null;
    lengthlimit := if opts.LengthLimit === infinity then 0 else opts.LengthLimit;
    if not instance(lengthlimit, ZZ) or lengthlimit < 0 then return null;
    m := presentation M;
    if not msolveApplicable m then return null;
    w := rawMsolveMinimalBetti(raw m, lengthlimit, allowableThreads, gbTrace);
    if w === null then return null;
    betti(unpackMsolveBetti w, Weights => opts.Weights))

-- msolveSetup() installs all of them; msolveSetup {gb, eliminate} a selection
msolveSetup = arg -> (
    install := if instance(arg, VisibleList) then toList arg else {arg};
    if #install == 0 then install = {
        gb, groebnerBasis, -* syz, *- mingens, trim,
        eliminate, kernel, saturate,
        ContainmentHooks, minimalBetti};
    printerr("installing msolve hooks for ", install);

    if not rawMsolvePresent() then printerr(
	"warning: this Macaulay2 was not built against the msolve library; ",
	"msolveSetup will route through the msolve executable instead");

    if member(gb, install) then
    gb Matrix := opts -> m -> (
        msolveGBHook(opts, m) ?? (M2DefaultGB opts)(m));

    if member(groebnerBasis, install) then
    groebnerBasis Matrix := opts -> m -> (
        generators msolveGBHook(opts, m) ?? (M2DefaultGBasis opts)(m));

    if member(syz, install) then
    addHook((syz, Matrix), Strategy => Msolve, msolveSyzygyHook);

    if member(eliminate, install) then
    eliminate(List, Ideal) := Ideal => (v, I) -> (
        R := ring I;
        if #v === 0 then return I;
        if any(v, x -> ring x =!= R)
        then error "expected a list of elements in the ring of the ideal";
        msolveEliminationGB(I, v) ?? M2DefaultEliminate(v, I));

    if member(kernel, install) then
    addHook((kernel, RingMap), Strategy => Msolve, (opts, f) -> (
	    (F, R) := (target f, source f);
	    if not isAffineRing R or not isAffineRing F
	    or coefficientRing R =!= coefficientRing F
	    or opts.?SubringLimit and opts.SubringLimit =!= infinity then return null;
	    -- the kernel is what survives eliminating the variables of the target from the graph ideal
	    g := generators graphIdeal f;
	    n1 := numgens F;
	    result := msolveElimGB(g, n1);
	    if result === null then return null;
	    (S1, G) := result;
	    mapback := map(R, ring g, map(R^1, R^n1, 0) | vars R);
	    ideal mapback substitute(G, vars ring g)));

    if member(saturate, install) then
    addHook((saturate, Ideal, RingElement), Strategy => Msolve,
	(opts, I, f) -> try ideal msolveSaturate(I, f));

    -- by default, Ideal == ZZ calls rawGBContains
    if member(ContainmentHooks, install) then
    addHook(ContainmentHooks, Strategy => Msolve,
        (f, g) -> try f % msolveGBHook g == 0);

    if member(mingens, install) then
    addHook((mingens, Module), Strategy => Msolve, (opts, M) -> (
            if M.?relations or not M.?generators
            or not msolveApplicable M.generators
            then return null;

            R := ring M;
            m := vars R;

            -- Find classes in M / m*M
            f := M.generators;
            -- msolve may decline the generators of m ** f even when it
            -- accepts f itself, in which case the hook defers to Macaulay2
            if (G := msolveGBHook(m ** f)) === null then return null;
            B := f % G;

            -- Select a kk-basis of the classes
            kk := coefficientRing R;
            -- TODO: this should be embarrasingly degree parallelized
            (N, C) := coefficients B;
            cols := columnRankProfile mutableMatrix lift(C, kk);
            N * C_cols));

    if member(trim, install) then
    addHook((trim, Module), Strategy => Msolve, (opts, M) -> (
            if M.?relations or not M.?generators
            or not msolveApplicable M.generators
            then return null;

            if (g := mingens M) === M.generators
            then M else image mingens M));

    if member(minimalBetti, install) then
    minimalBetti Module := opts -> M -> (
        msolveMinimalBettiHook(opts, M) ?? (M2DefaultMinimalBetti opts)(M));

    ///
    -- TODO: what other rawGB... compiled functions can use msolve?
    -- TODO: try out pushforward using msolve?

    if member(ReduceHooks, install) then
    if member(length, install) then
    if member(prune, install) then
    if member(pushForward, install) then
    if member(quotient, install) then
    if member(quotientRemainder, install) then
    if member(remainder, install) then
    if member(res, install) then
    ///;
    )

--------------------------------------------------------------------------------
-- Rational interval type, constructors, and basic methods
--------------------------------------------------------------------------------

QQi = new Ring of List -- TODO: array looks better, but List implements arithmetic by default!
QQi.synonym = "rational interval"

ring      QQi := x -> QQi
precision QQi := x -> infinity

QQinterval = method(TypicalValue => QQi)
QQinterval VisibleList := bounds -> (
    if #bounds == 2 then QQinterval(bounds#0, bounds#1)
    else error "expected a lower bound and upper bound")
QQinterval Number                        := midpt  -> QQinterval(midpt/1, midpt/1)
QQinterval InexactNumber                 := midpt  -> QQinterval lift(midpt, QQ)
QQinterval(InexactNumber, InexactNumber) := (L, R) -> QQinterval(lift(L, QQ), lift(R, QQ))
QQinterval(Number,        Number)        := (L, R) -> new QQi from [L/1, R/1]

-- TODO: these are compiled functions, make them methods and define for QQi
left'     = first
right'    = last
midpoint' = int -> sum int / 2
diameter QQi := x -> x#1 - x#0

interval QQi := opts -> x -> interval(x#0, x#1, opts)

QQi == Number :=
Number == QQi := (x, y) -> QQinterval x == QQinterval y

-- ZZ, QQ, RR, RRi
promote(Number, QQi) := (n, QQi) -> QQinterval n

-- ZZ, QQ
lift(QQi, Number) := o -> (x, R) -> (
    if diameter x == 0 then lift(midpoint' x, R)
    else if o.Verify then error "lift: interval has positive diameter")

--------------------------------------------------------------------------------

msolveDefaultPrecision = 32 -- alternative: defaultPrecision

msolveRealSolutions = method(TypicalValue => List, Options => msolveDefaultOptions)
msolveRealSolutions Ideal              := opt ->  I0     -> msolveRealSolutions(I0, QQi, opt)
msolveRealSolutions(Ideal, RingFamily) := opt -> (I0, F) -> msolveRealSolutions(I0, F_msolveDefaultPrecision, opt)
msolveRealSolutions(Ideal, Ring)       := opt -> (I0, F) -> (
    if not any({QQ, QQi, RR_*, RRi_*}, F' -> ancestor(F', F))
    then error "msolveRealSolutions: expected target field to be rationals, reals, or a rational or real interval field";
    (S, K, I) := toMsolveRing I0;
    -- if precision is not specified, we want to use msolve's default precision
    prec := if precision F === infinity then msolveDefaultPrecision else precision F;
    mOut := msolve(S, K, I_*, "-p " | prec, opt);
    -- format: [dim, [numlists, [ solution boxes ]]] when zero-dimensional, otherwise [1?, numgens, -1, []]
    mSeq := toSequence value readMsolveList get mOut;
    d := mSeq#0;
    if d =!= 0 then error "msolveRealSolutions: expected zero dimensional system of equations";
    solsp := mSeq#1;
    if solsp_0 > 1 then (
	printerr "msolveRealSolutions: unexpected msolve output, returning full output"; return {d, solsp});
    prec  = max(defaultPrecision, prec); -- we want the output precision to be at least defaultPrecision
    sols := apply(last solsp, sol -> apply(sol, QQinterval));
    if ancestor(QQi,   F) then sols else
    if ancestor(QQ,    F) then apply(sols, sol -> apply(sol, midpoint'))                    else
    if ancestor(RR_*,  F) then apply(sols, sol -> apply(sol, midpoint') * numeric(prec, 1)) else
    if ancestor(RRi_*, F) then apply(sols, sol -> apply(sol, range -> interval(range, Precision => prec))))

msolveRUR = method(TypicalValue => List, Options => msolveDefaultOptions)
msolveRUR Ideal := opt -> I0 ->(
    S0 := ring I0;
    (S, K, I) := toMsolveRing I0;
    mOut := msolve(S, K, I_*, "-P 2", opt);
    -- format: [dim, [char, nvars, deg, vars, form, [1, [lw, lwp, param]]]]:
    solsp := value readMsolveList get mOut;
    if first solsp != 0 then error "msolveRUR: expected zero dimensional input ideal";
    lc:=(solsp_1)_4;
    l:=sum(numgens S0,i->lc_i*S0_i);
    RUR:= new MutableHashTable;
    T:= getSymbol("T");
    S2 := K(monoid[T]);
    T=first gens S2;
    RUR#"T"=l;
    RUR#"var"=T;
    RUR#"degree"=(solsp_1)_2;
    para:= ((solsp_1)_5)_1;
    sVarsMsolve:=(solsp_1)_3;
    W:=sum((para_0)_0+1,i->(T)^i*(((para_0)_1)_i));
    RUR#"findRootsUniPoly"=W;
    RUR#"denominator"=diff(T,W);
    vs:=last para;
    msolveVarOrder:=for sg in gens(S) list position(sVarsMsolve,i->i==toString(sg));
    tempNumerator:=(append(for f in vs list sum((f_0)_0+1,i->T^i*((f_0)_1)_i),-T*diff(T,W)));
    if (numgens(S)==#(tempNumerator)) then(
	RUR#"numerator"=tempNumerator_msolveVarOrder;
	)
    else(
	RUR#"numerator"=tempNumerator_(append(msolveVarOrder, numgens(S)));
	);
    return new HashTable from RUR;
    );

load "./Msolve/tests.m2"

beginDocumentation()

doc ///
Node 
     Key
     	  Msolve
     Headline
	  Macaulay2 interface for msolve; computes real solutions and Groebner basis, etc.
     Description
     	  Text
              This package provides a Macaulay2 interface for the
	      msolve library [1] developed by
              Jérémy Berthomieu, Christian Eder, and Mohab Safey El
              Din.
	      
	      The package has functions to compute Groebner basis, in
	      @TO GRevLex@ order only, for ideals with rational or finite
	      field coefficients. Finite field characteristics must be
	      less than $2^{31}$. There are also functions to
	      compute elimination ideals, for ideals with rational or
	      finite field coefficients.
	      
	      The @TO2 {"Saturation::saturate", "saturation"}@ of an ideal by a single polynomial may be
	      computed for ideals with finite field coefficients, again
	      with characteristic less than $2^{31}$.
	      
	      For zero dimensional polynomial ideals, with integer or
	      rational coefficients, there are functions to compute all
	      real solutions, and to compute a rational univariate
	      representation of all (complex) solutions.
	      
	      The M2 interface assumes that the binary executable is
	      named "msolve" is on the executable path.

	      For all functions the option @TT "Verbosity"@ can be used.
	      It has levels 0, 1, 2. The default is 0.

	      Msolve supports parallel computations. The option @TT "Threads"@ is used to set this.
	      The default value is allowableThreads, but this can be set manually by the user when 
	      calling a function. E.g. for an ideal I:
	  Example
	      R = QQ[x,y,z]
	      I = ideal(x, y, z)
	      msolveGB(I, Verbosity => 2, Threads => 6) -* no-capture-flag *-
    References
      [1] The msolve library: @HREF {"https://msolve.lip6.fr"}@;

Node 
    Key
    	msolveGB
       (msolveGB, Ideal)
       [msolveGB, Threads]
       [msolveGB, Verbosity]
       [msolveGB, DegreeLimit]
    Headline
	compute generators of a Groebner basis in GRevLex order
    Usage
    	msolveGB(I)
    Inputs
    	I:Ideal
	    in a polynomial ring with @TO GRevLex@ order and coefficients over @TO QQ@ or
	    @TO2 {"finite fields", TT "ZZ/p"}@ in characteristic less than $2^{31}$
	Threads => ZZ -- number of processor threads to use
	Verbosity => ZZ -- level of verbosity between 0, 1, and 2
	DegreeLimit => List -- stop after this degree; the default @TT "{}"@ computes the whole basis
    Outputs
        GB:Matrix
	    whose columns form a Groebner basis for the input ideal I, in the GRevLex order
    Description 
        Text
	    This functions uses the F4 implementation in the msolve package to compute a Groebner basis,
	    in GRevLex order, of a polynomial ideal with either rational coefficients or finite field
	    coefficients with characteristic less than $2^{31}$. If the input ideal is a polynomial ring
	    with monomial order other than GRevLex a GRevLex basis is returned (and no warning is given).
	    The input ideal may also be given in a ring with integer coefficients, in this case a Groebner
	    basis for the given ideal over the rationals  will be computed, denominators will be cleared,
	    and the output will be a Groebner basis over the rationals in GRevLex order with integer coefficients.
    	Text
	    First an example over a finite field
	Example
	    R=ZZ/1073741827[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    gB=msolveGB I
	    lT=monomialIdeal leadTerm gB
	    degree lT
	    dim lT	    
	Text
	    Now the same example over the rationals. 
	Example 
	    R=QQ[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    gB=msolveGB I
	    (ideal gB)== ideal(groebnerBasis I)
	    lT=monomialIdeal leadTerm gB
	    degree lT
	    dim lT
	Text
	    A @TO DegreeLimit@ stops msolve's degree by degree schedule once every
	    S-pair of that degree has been handled.  The degree may be a
	    multidegree, in which case it is collapsed to a single integer against
	    the ring's heft vector, exactly as @TO gb@ does before Macaulay2's own
	    engine sees one; everything of no greater heft is computed as well.
	Example
	    R = ZZ/32003[x,y,z]
	    I = ideal(x^2+y*z, y^2+x*z, z^2+x*y)
	    numcols msolveGB(I, DegreeLimit => {2})
	    numcols msolveGB I
	Text
	    What comes back is a Groebner basis of nothing in particular: it
	    generates the ideal only in degrees up to the limit.  Nor can it be
	    continued -- msolve keeps no state between calls, so asking again for a
	    larger limit recomputes from the input.  For that reason a limited
	    request returns the matrix rather than a @TO GroebnerBasis@, and
	    @TO msolveSetup@ does not route @TO gb@'s own @TO DegreeLimit@ here.
Node
    Key
        msolveSyzygy
       (msolveSyzygy, Matrix)
       [msolveSyzygy, Threads]
       [msolveSyzygy, Verbosity]
       [msolveSyzygy, DegreeLimit]
       [msolveSyzygy, SyzygyLimit]
       [msolveSyzygy, SyzygyRows]
    Headline
        compute syzygies of the columns of a matrix using msolve
    Usage
        msolveSyzygy M
    Inputs
        M:Matrix
            over a GRevLex polynomial ring over a prime field of characteristic less than $2^{31}$
        Threads => ZZ -- number of processor threads to use
        Verbosity => ZZ -- level of verbosity between 0, 1, and 2
        DegreeLimit => List -- stop after this degree; the default @TT "{}"@ computes everything
        SyzygyLimit => ZZ -- stop after this many syzygies; the default @TO infinity@ computes all of them
        SyzygyRows => ZZ -- keep only this many rows; the default @TO infinity@ keeps all of them
    Outputs
        S:Matrix
            whose columns generate the syzygy module of the columns of M
    Description
        Text
            This function uses msolve's graph-module computation.  The result
            is a Groebner basis of the syzygy module and need not be a minimal
            generating set.
        Example
            R = ZZ/32003[x,y,z,w]
            M = matrix {{x,y,z}, {y,z,w}}
            S = msolveSyzygy M
            M * S == 0
            image S == image syz M
        Text
            @TO SyzygyLimit@ is a genuine early stop rather than a cap on the
            output: under the position over term order a basis element is a
            syzygy exactly when its lead term is, so msolve counts them as it
            goes and quits.  A degree is always finished, so more relations than
            asked for can turn up at once, and the extras are discarded.
        Example
            N = matrix {{x,y,z}}
            numcols msolveSyzygy(N, SyzygyLimit => 2)
            numcols msolveSyzygy N
        Text
            @TO SyzygyRows@ keeps only the first rows of the syzygy matrix and
            drops the columns that become zero.  The result is a submatrix of
            the syzygy matrix, so it is @EM "not"@ annihilated by M; it is what a
            caller wants when the columns of M are generators followed by
            relations and only the coefficients on the generators matter.
        Example
            T = msolveSyzygy(N, SyzygyRows => 1)
            N * T == 0
Node
    Key
        msolveResolution
       (msolveResolution, Matrix)
       (msolveResolution, Ideal)
       (msolveResolution, Module)
       [msolveResolution, Threads]
       [msolveResolution, Verbosity]
       [msolveResolution, LengthLimit]
        MsolveResolution
        MsolveDifferential
       (symbol _, MsolveResolution, ZZ)
       (symbol _, MsolveDifferential, ZZ)
       (length, MsolveResolution)
       (ring, MsolveResolution)
       (net, MsolveResolution)
       (complex, MsolveResolution)
    Headline
        a free resolution computed by msolve, one differential at a time
    Usage
        C = msolveResolution M
    Inputs
        M:{Matrix,Ideal,Module}
            over a GRevLex polynomial ring over a prime field of characteristic
            less than $2^{31}$, homogeneous, and singly graded with each
            variable's degree equal to its weight
        Threads => ZZ -- number of processor threads to use
        Verbosity => ZZ -- level of verbosity between 0, 1, and 2
        LengthLimit => ZZ -- truncate the resolution at this level
    Outputs
        C:MsolveResolution
    Description
        Text
            Unlike the other functions here, this one returns a live
            computation rather than an answer.  Building it runs the module
            Groebner basis and the whole Schreyer frame.  The frame is
            combinatorial -- no field arithmetic happens past the Groebner
            basis -- so from that point on every free module @TT "C_i"@ in the
            resolution, its rank and the degrees of its generators, is free to
            ask for.
        Example
            R = ZZ/32003[x,y,z,w]
            I = minors_2 matrix {{x,y,z}, {y,z,w}}
            C = msolveResolution I
            length C
            C_2
            degrees C_2
        Text
            Asking for a differential is what makes msolve reduce, and then
            only up to the level asked for: the block at level $i$ in degree
            $d$ reduces against level $i-1$ in degrees at most $d$, so $d_i$
            needs $d_2, \dots, d_{i-1}$ and nothing above.  Levels already
            computed are remembered, so asking twice costs nothing and asking
            out of order costs no more than asking in order.
        Example
            C.dd_2
            C.dd_1 * C.dd_2 == 0
            image C.dd_1 == module I
        Text
            The complex is the @EM "nonminimal"@ one, as with
            @TT "res(..., Strategy => Nonminimal)"@: @TT "C_1"@ is the Groebner
            basis of the image of the input, not the input columns, since
            msolve keeps no change of basis between the two.  Its ranks
            therefore depend on the Groebner basis, hence on the module order,
            and need not agree with Macaulay2's own nonminimal resolution of
            the same module.  Only the minimal Betti numbers are an invariant.
        Text
            @TT "LengthLimit"@ truncates the frame.  Note there is no level
            that is always enough: the frame is nonminimal, so Hilbert's syzygy
            theorem does not bound it, and it really can run past the number of
            variables.
        Example
            S = ZZ/32003[x,y,z]
            msolveResolution ideal(z, y^2, x^2*y, x^3)
            msolveResolution(ideal(z, y^2, x^2*y, x^3), LengthLimit => 2)
        Text
            @TT "complex C"@ materializes every differential and hands back an
            ordinary @TO Complex@, which is the opposite of what the live handle
            is for and so is never done implicitly.  What it buys is
            @TO (minimize, Complex)@: the nonminimal resolution msolve computes
            becomes the minimal one, resolving the cokernel of the input.
        Example
            use R
            D = complex C
            betti minimize D
            betti minimize D == minimalBetti I
        Text
            That is the expensive way to get a minimal Betti table, and it is
            here to check the cheap one rather than to compete with it: rank
            extraction reads the same numbers off the ranks of the scalar parts
            of the differentials, never building a minimal complex at all.
            With a @TT "LengthLimit"@ the top level is the exception -- its
            minimal rank is only correct once the level above it has been seen,
            so minimizing a truncated complex can overstate the last one.
Node
    Key
    	msolveLeadMonomials
       (msolveLeadMonomials, Ideal)
       [msolveLeadMonomials, Threads]
       [msolveLeadMonomials, Verbosity]
    Headline
	compute the leading monomials of a Groebner basis in GRevLex order
    Usage
    	msolveLeadMonomials(I)
    Inputs
    	I:Ideal
	    in a polynomial ring with @TO GRevLex@ order and coefficients over @TO QQ@ or
	    @TO2 {"finite fields", TT "ZZ/p"}@ in characteristic less than $2^{31}$
	Threads => ZZ -- number of processor threads to use
	Verbosity => ZZ -- level of verbosity between 0, 1, and 2
    Outputs
        GB:Matrix
	    whose columns are the leading monomials (of a Groebner basis for) the input ideal I, in the GRevLex order
    Description 
        Text
	    This functions uses the F4 implementation in the msolve package to compute leading
	    monomials via a Groebner basis, in GRevLex order, of a polynomial ideal with either
	    rational coefficients or finite field coefficients with characteristic less than $2^{31}$.
	    If the input ideal is a polynomial ring with monomial order other than GRevLex a GRevLex
	    basis is returned (and no warning is given). The input ideal may also be given in a ring
	    with integer coefficients, in this case a Groebner basis for the given ideal over the
	    rationals  will be computed, denominators will be cleared, and the output will be a
	    Groebner basis over the rationals in GRevLex order with integer coefficients.
    	Text
	    First an example over a finite field
	Example
	    R=ZZ/1073741827[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    lm=monomialIdeal msolveLeadMonomials I
	    degree lm
	    dim lm	    
	Text
	    Now the same example over the rationals; note over the rationals msolve first
	    computes over a finite field and when only the leading monomials are asked for
	    the correct leading monomials will be returned but the full Groebner basis over
	    @TO QQ@ will not be computed. Hence if only degree and dimension are desired
	    this command will often be faster that the Groebner basis command.
	Example 
	    R=QQ[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    lm=monomialIdeal msolveLeadMonomials I
	    lt=monomialIdeal leadTerm groebnerBasis I
	    lm==lt
	    degree lm
	    dim lm

Node
    Key
        QQi
    Headline
        the class of all rational intervals
    Description
        Text
            This class is similar to the class of @TO2(RRi, "real intervals")@,
	    except that the boundaries are arbitrary precision rational numbers.
    Caveat
        Currently this class is not implemented in the interpreter,
	which means rings or matrices over rational intervals are not supported,
	and many functionalities of @TO RRi@ are not yet available.

Node 
    Key
    	msolveRealSolutions
       (msolveRealSolutions, Ideal)
       (msolveRealSolutions, Ideal, Ring)
       (msolveRealSolutions, Ideal, RingFamily)
       [msolveRealSolutions, Threads]
       [msolveRealSolutions, Verbosity]
    Headline
	compute all real solutions to a zero dimensional system using symbolic methods
    Usage
    	msolveRealSolutions(I)
	msolveRealSolutions(I, K)
    Inputs
    	I:Ideal
	    which is zero dimensional, in a polynomial ring with coefficients over @TO QQ@
	K:{Ring,RingFamily}
	    the field to find answers in, which must be one of
	    @TO QQi@ (default), @TO QQ@, @TO RR@, or @TO RRi@ (possibly with specified precision)
	Threads => ZZ -- number of processor threads to use
	Verbosity => ZZ -- level of verbosity between 0, 1, and 2
    Outputs
        :List
	    of lists; each entry in the list consists of a list representing
	    the coordinates of a solution. By default each solution coordinate value is
	    also represented by a @TO2(QQi, "rational interval")@ consisting of a two element
	    list of rational numbers, @TT "{a, b}"@, this means that that coordinate of the
	    solution has a value greater than or equal to @TT "a"@ and less than or equal to @TT "b"@.
	    This interval is computed symbolically and its correctness is guaranteed by exact methods.
    Description 
        Text
	    This functions uses the msolve package to compute the real solutions to a zero
	    dimensional polynomial ideal with either integer or rational coefficients.
	Text
	    The second input is optional, and indicates the alternative ways to provide output
	    either using an exact rational interval @TO QQi@, a real interval @TO RRi@,
	    or by taking a rational or real approximation of the midpoint of the intervals.
	Example
	    R = QQ[x,y]
	    I = ideal {(x-1)*x, y^2-5}
	    rationalIntervalSols = msolveRealSolutions I
	    rationalApproxSols = msolveRealSolutions(I, QQ)
	    floatIntervalSols = msolveRealSolutions(I, RRi)
	    floatIntervalSols = msolveRealSolutions(I, RRi_10)
	    floatApproxSols = msolveRealSolutions(I, RR)
	    floatApproxSols = msolveRealSolutions(I, RR_10)
	Text
	    Note in cases where solutions have multiplicity this is not reflected in the output.
	    While the solver does not return multiplicities,
	    it reliably outputs the verified isolating intervals for multiple solutions.  
	Example 
	    I = ideal {(x-1)*x^3, (y^2-5)^2}
	    floatApproxSols = msolveRealSolutions(I, RRi)
Node 
    Key
    	msolveRUR
       (msolveRUR, Ideal)
       [msolveRUR, Threads]
       [msolveRUR, Verbosity]
    Headline
	compute the rational univariate representation using symbolic methods
    Usage
    	msolveRUR(I)
    Inputs
    	I:Ideal
	    which is zero dimensional, in a polynomial ring with coefficients over @TO QQ@
	Threads => ZZ -- number of processor threads to use
	Verbosity => ZZ -- level of verbosity between 0, 1, and 2
    Outputs
        RUR:HashTable
	    with 6 keys giving the rational univariate representation of I
    Description 
        Text
	    This functions uses the msolve package to compute the rational univariate representation
	    (RUR) of a zero dimensional polynomial ideal with either integer or rational coefficients.
	    
	    The RUR gives a parametrization for all complex solutions to the input system.
	    For a complete definition of the RUR see the paper: Rouillier, Fabrice (1999).
	    "Solving Zero-Dimensional Systems Through the Rational Univariate Representation".
	    Appl. Algebra Eng. Commun. Comput. 9 (9): 433–461.
	    
	    If I is a zero dimensional ideal in QQ[x_1..x_n] then the RUR is given by:
	    
	    (x_1,..,x_n)={ (-v_1(T)/w'(T), .. , -v_n(T)/w'(T)) | w(T)=0}
	    
	    The output is a hash table with 6 keys. 
	    
	    The key "degree" is the number of solutions to I, counted with multiplicity. 
	    
	    The key "findRootsUniPoly" gives the polynomial w(T) above.
	    
	    The key "denominator"  gives the polynomial w'(T), which is the derivative of w(T) and
	    is the denominator of each coordinate above.
	    
	    The key "numerator" gives a list {v_1(T), .. , v_n(T)} of length n above, with n the
	    number of variables, where the polynomial v_i(T) gives the numerator of the ith coordinate.
	    
	    The key "var" gives the variable name in the univariate polynomial ring; by default this is: "T".
	    
	    The key "T" gives the linear relation between the variables of the ring of I and
	    the single variable, which is denoted T above.
	    
    	Text
	    A simple example, where the input ideal is zero dimensional and radical.
	Example
	    R = QQ[x_1..x_3]
	    f = (x_1-1)
	    g = (x_2-2)
	    h = (x_3^2-9)
	    I = ideal (f,g,h)
	    decompose I
	    rur=msolveRUR(I)
	    factor rur#"findRootsUniPoly"
	    sols=-1*(rur#"numerator")	    
	    denom= rur#"denominator"
	    (for s in sols list sub(s,T=>3))/sub(denom,T=>3)
	    (for s in sols list sub(s,T=>-3))/sub(denom,T=>-3)	    
	Text
	    In cases where the input ideal has dimension greater than zero an error will be returned.     

Node 
    Key
    	msolveSaturate
       (msolveSaturate, Ideal,RingElement)
       [msolveSaturate, Threads]
       [msolveSaturate, Verbosity]
    Headline
	compute a Groebner basis for the saturation of an ideal by a single polynomial in GRevLex order
    Usage
    	msolveSaturate(I)
    Inputs
    	I:Ideal
	    in a polynomial ring with @TO GRevLex@ order and coefficients over
	    @TO2 {"finite fields", TT "ZZ/p"}@ in characteristic more than $2^16$ and less than $2^{31}$
	f:RingElement
	    a polynomial in the same ring as I.    
	Threads => ZZ -- number of processor threads to use
	Verbosity => ZZ -- level of verbosity between 0, 1, and 2
    Outputs
        GB:Matrix
	    whose columns form a Groebner basis for the ideal $I:f^\infty$, in the GRevLex order
    Description 
        Text
	    This functions uses the F4SAT algorithm implemented in the msolve library to compute a
	    Groebner basis, in GRevLex order, of $I:f^\infty$, that is of the saturation of the
	    ideal $I$ by the principal ideal generated by the polynomial $f$.
	Example
	    R=ZZ/1073741827[z_1..z_3]
	    I=ideal(z_1*(z_2^2-17*z_1-z_3^3),z_1*z_2)
	    satMsolve=ideal msolveSaturate(I,z_1)
	    satM2=saturate(I,z_1)
	Text
	    Note that the ring must be a polynomial ring over a finite field.
    Caveat
	Currently the F4SAT algorithm is only implemented over prime fields in characteristic between $2^{16}$ and $2^{31}$.

Node 
    Key
    	msolveEliminate
	(msolveEliminate, Ideal,List)
	(msolveEliminate, List,Ideal)
	(msolveEliminate, Ideal,RingElement)
	(msolveEliminate, RingElement,Ideal)
       [msolveEliminate, Threads]
       [msolveEliminate, Verbosity]
    Headline
	compute the elimination ideal of a given ideal
    Usage
    	msolveEliminate(I,elimVars)
    Inputs
    	I:Ideal
	    in a polynomial ring with @TO GRevLex@ order and coefficients over
	    @TO2 {"finite fields", TT "ZZ/p"}@ in characteristic less than $2^{31}$
	elimVars:List
	    of variables in the same ring as @TT "I"@, these variables will be eliminated
	Threads => ZZ -- number of processor threads to use
	Verbosity => ZZ -- level of verbosity between 0, 1, and 2
    Outputs
        GB:Matrix
	    whose columns are generators for the elimination ideal
    Description 
        Text
            This function takes as input a polynomial ideal and
            computes the elimination ideal given by eliminating the
            variables specified in the inputted list.

            The behavior is very different over finite (prime) fields, and the rationals.
	    Over @TO QQ@, the subideal over a smaller set of variables eliminating the given ones is returned.
            Over a finite field, the Groebner basis in the product order eliminating the given block of variables
	    is returned (warning: this is a copy of the ring with potentially permuted variables).
        Text
	    First an example over the rationals. 
	Example
	    R = QQ[x,a,b,c,d]
	    f = 7*x^2+a*x+b
	    g = 2*x^2+c*x+d
	    M2elim=eliminate(x,ideal(f,g))
	    Msolveelim=msolveEliminate(x,ideal(f,g))
	    sub(M2elim,ring Msolveelim)==Msolveelim    
	Text
            We can also work over a finite field. Here we get the full
            Groebner basis in the permuted variables with a block order.
        Example 
	    R = ZZ/1073741827[x,y,a,b,c,d]
	    f = c*x^2+a*x+b*y^2+d*x*y+y
	    g = diff(x,f)
	    h = diff(y,f)
	    M2elim = eliminate({x,y}, ideal(f,g,h))
	    Msolveelim = msolveEliminate({x,y}, ideal(f,g,h))
	    M2elim_0 == sub(Msolveelim_0, R)
///	      

-*
    uninstallPackage "Msolve"
    restart
    needsPackage("Msolve")
    installPackage "Msolve"
    restart
    check "Msolve"
*-
