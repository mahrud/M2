--		Copyright 1996-2002 by Daniel R. Grayson

needs "enginering.m2"
needs "quotring.m2"
needs "orderedmonoidrings.m2"

protect FieldExtension

-----------------------------------------------------------------------------
-- Local variables
-----------------------------------------------------------------------------

-- Flint's implementation of Zech logarithm representation
-- for Galois fields assumes p^d < 2^FLINT_BITS
-- c.f. http://flintlib.org/sphinx/fq_zech.html#c.fq_zech_ctx_init
flintbits := version#"flint size"

-- Labels for algorithms implemented in the engine
rawGFtypes := new HashTable from {
    "Old"      => rawGaloisField,
    "New"      => rawARingGaloisFieldM2,
    "Flint"    => rawARingGaloisFieldFlintZech,
    "FlintBig" => rawARingGaloisFieldFlintBig, -- does not need a primitive element
    }

-----------------------------------------------------------------------------
-- Local utilities
-----------------------------------------------------------------------------

packGF := (R, p, d, f) -> (
    kk := R/f;
    kk.char = p;
    kk.degree = d;
    kk.order = p^d;
    kk)

unpack := S -> (			  -- a quotient ring
    if not instance(S, QuotientRing)
    then error("GF: expected a quotient ring");
    R := ultimate(ambient, S);			  -- original poly ring
    if R === ZZ
    then if isField S then ( R = S[local t]; S = R/t; )
    else error("GF: expected the base ring to be a field");
    --
    if not numgens R === 1
    or not instance(R, PolynomialRing)
    then error("GF: expected a polynomial ring in one variable");
    -- why?
    if degreeLength R =!= 1
    then error("GF: expected a polynomial ring with standard grading");
    --
    p := char R;
    if not isPrime p
    or not coefficientRing R === ZZ/p
    then error("GF: expected a polynomial ring over a finite prime field");
    --
    I := ideal S;			  -- a Matrix, sigh
    if not numgens I === 1
    or not isPrime(f := I_0)
    then error("GF: expected a quotient ring by an irreducible polynomial");
    --
    d := first degree f;
    (R, p, d, f))

getVariable := (a, defaultVar) -> (
    if instance(a, Nothing) then defaultVar else
    if instance(a, String) then getSymbol a else baseName a)

ZZpA = (p, a) -> (ZZ/p) (monoid [ getVariable(a, getSymbol "a") ])
ZZpA = memoize ZZpA -- cached globally, a common base ring for the extension fields

-----------------------------------------------------------------------------
-- GaloisField type declarations, constructors, and basic methods
-----------------------------------------------------------------------------

GaloisField = new Type of EngineRing
GaloisField.synonym = "Galois field"

-- constructors
GF = method (
    Options => {
	PrimitiveElement => FindOne, -- TODO: deprecate FindOne
	SizeLimit => 2^16,
	Strategy  => null, -- see hooks(GF, ZZ, ZZ) for strategies
	Variable  => null,
	}
    )

GF ZZ := GaloisField => opts -> q -> (
    pd := factor q;
    if #pd =!= 1 or pd#0#0 === -1
    then error "GF: expected a power of a prime";
    GF(pd#0#0, pd#0#1, opts))

GF(ZZ, ZZ) := GaloisField => opts -> (p, d) -> (
    if not isPrime p then error "GF: expected a prime number as base";
    if d < 0         then error "GF: expected non-negative exponent";
    (a, type) := (opts. Variable, opts.Strategy);
    -- return a ZZ/p field
    if d == 1 and member(type, keys rawZZptypes) then return ZZp(p, Strategy => type);
    -- find an irreducible polynomial
    -- note: f is an ideal so that S/f is cached
    S := ZZpA(p, a); -- ZZ/p[a]
    -- TODO: use https://github.com/Macaulay2/M2/issues/1596 when it is implemented
    f := if S#?(FieldExtension, d) then S#(FieldExtension, d)
    else S#(FieldExtension, d) = ideal runHooks((GF, ZZ, ZZ), Strategy => type, (p, d, S_0));
    -- find a primitive generator
    g := findPrimitive packGF(S, p, d, f);
    GF(ring g, opts ++ { PrimitiveElement => g, Variable => S_0 }))

-- TODO: use https://github.com/Macaulay2/M2/issues/1596 when it is implemented
GF Ring := GaloisField => opts -> S -> if S#?(symbol GF, opts) then return S#(symbol GF, opts) else S#(symbol GF, opts) = (
    (R, p, d, f) := unpack S;
    (g, type, a) := (
	opts.PrimitiveElement,
	opts.Strategy,
	opts.Variable);
    --
    if g === FindOne then g = findPrimitive S else (
	if not ring g === S  then error "GF: expected primitive element in the right ring";
	if not isPrimitive g then error "GF: expected ring element to be primitive");
    a = getVariable(a, S.generatorSymbols#0);
    --
    if type === null then type = "Flint";
    if not rawGFtypes#?type then error("GF: unrecognized type ", toString type, " of Galois field requested;", newline,
	"  available types are: ", demark(", ", keys rawGFtypes));
    -- TODO: why does it make sense to forget the primitive generator if we already have it?
    if opts.SizeLimit < p^d or g =!= S_0 then (
	type = "FlintBig";
	g = S_0; -- Possibly NOT the primitive element in this case!!  We don't need it, and we don't want to compute it yet.
	-- Note: we used to have Galois fields always encoded by powers of the primitive element.  Ring maps would use this
	-- (in ringmap.m2) to help tell the engine where the primitive element goes.
	-- But: for FlintBig, this isn't being used.  We should perhaps consider changing the code in ringmap.m2.
	-- For now, setting g will do.
	);
    --
    rawF := rawGFtypes#type raw g;
    F := new GaloisField from rawF;
    F.PrimitiveElement = g; -- Note: the primitive element is in S
    F.baseRings = append(S.baseRings, S);
    F.degreeLength = 0;
    -- TODO: what are these used for?
    F.rawGaloisField = true;
    F.isBasic = true;
    F.promoteDegree = makepromoter 0;			    -- do this before commonEngineRingInitializations
    F.liftDegree = makepromoter degreeLength S;	    -- do this before commonEngineRingInitializations
    commonEngineRingInitializations F;
    --
    F.isCommutative = true;
    F.char = p;
    F.degree = d;
    F.order = p^d;
    F.frac = F;
    F.generators = apply(generators S, m -> promote(m,F)); -- this will be wrong if S is a tower
    --
    if S.?generators           then F.generators           = apply(S.generators, r -> promote(r,F));
    if S.?generatorSymbols     then F.generatorSymbols     = S.generatorSymbols;
    if S.?generatorExpressions then F.generatorExpressions = S.generatorExpressions;
    if S.?indexSymbols then F.indexSymbols = applyValues(S.indexSymbols, r -> promote(r,F));
    if S.?indexStrings then F.indexStrings = applyValues(S.indexStrings, r -> promote(r,F));
    --
    baseName   F := y -> if F_0 == y then a else error "expected a generator";
    expression F := t -> expression lift(t, S);
    --
    F / F := (x,y) -> if y == 0 then error "division by zero" else x // y;
    F % F := (x,y) -> if y == 0 then x else 0_F;
    use F)

-- basic methods
dim GaloisField := R -> 0

precision GaloisField := F -> infinity

generators GaloisField := opts -> F -> (
    if opts.CoefficientRing === null then F.generators
    else if opts.CoefficientRing === F then {}
    else apply(generators(baseRing F, opts), m -> promote(m, F)))

ambient GaloisField := Ring => baseRing

-- property checks
isField Ring := R -> R.?isField and R.isField

isAffineRing = method(TypicalValue => Boolean)
isAffineRing Ring := isField
isAffineRing PolynomialRing := R -> isCommutative R and not (options R).Inverses and isAffineRing coefficientRing R
isAffineRing QuotientRing := R -> isField R or isAffineRing ambient R

-- printing
expression GaloisField := F -> if hasAttribute(F, ReverseDictionary) then expression getAttribute(F, ReverseDictionary) else (expression GF) (expression F.order)
describe   GaloisField := F -> Describe (expression GF) (expression F.order)

net              GaloisField :=      net @@ expression
toString         GaloisField := toString @@ expression
toExternalString GaloisField := toString @@ describe

-- miscellaneous
random GaloisField := RingElement => opts -> F -> if (e := random F.order) > 0 then F_0^e else 0_F

isPrimitive = method(TypicalValue => Boolean)
isPrimitive RingElement := g -> isPrimitive(ring g, g)
isPrimitive(GaloisField,  RingElement) := (A, g) -> g != 0 and (
    all(first \ toList factor(qm1 := A.order - 1), p -> 1 != g ^ (qm1 // p)))
isPrimitive(QuotientRing, RingElement) := (R, g) -> g != 0 and (
    (S, p, d, f) := unpack R; qm1 := p^d-1;
    all(first \ toList factor qm1, v -> 1 != g ^ (qm1 // v)))

-- TODO: export this
findPrimitive = S -> (
    -- S is a polynomial ring of the form (ZZ/p)[x]/(f)
    (p, d, a) := (char S, degree S, S_0);
    while not isPrimitive a
    do (a = sum(d, i -> random p * a^i)); a)

-----------------------------------------------------------------------------
-- Strategies for finding a primitive element
-----------------------------------------------------------------------------

addHook((GF, ZZ, ZZ), Strategy => Default,
    (p, d, a) -> (
	f := a^d;
	while not isPrime f
	do (f = f + sum(d, i -> random p * a^i)); f))

addHook((GF, ZZ, ZZ), Strategy => "Flint",
    (p, d, a) -> if p^d < 2^(flintbits/2) then ( -- FIXME: use up to 2^FLINT_BITS
	-- TODO: is there a benefit to caching this somewhere?
	rawkk := rawARingGaloisField(p, d);
	c := rawARingGFPolynomial rawkk;
	sum(#c, i -> c#i * a^i)))

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/m2 "
-- End:
