-- Correctness pass.  Every method for a case runs in one process against input
-- built from a single ring: ZZ/p[w,x,y,z] evaluated twice gives two distinct
-- rings, and modules over them do not compare.  Each method still gets its own
-- module object, so a result cached by one does not answer for the next.
--
-- For elimination the two methods return their answers in different rings
-- (msolve's block order ring has the eliminated variables first), so the
-- comparison maps msolve's answer across by variable name.

needsPackage "Msolve"
debug Msolve
load "cases.m2"

bmFamily = getenv "FAMILY"
bmCaseNm = getenv "CASE"

bmCopy = M -> if M.?relations then subquotient(generators M, relations M) else image generators M

bmByName = (R, S) -> (
    nms := apply(gens R, toString);
    map(R, S, apply(gens S, v -> (
		i := position(nms, n -> n == toString v);
		if i === null then error("no variable named " | toString v);
		R_i))))

bmSay = (agree, detail) -> (
    << "VERIFY," << bmFamily << "," << bmCaseNm << ","
    << (if agree then "agree" else "DIFFER") << "," << detail << endl << flush)

bmCheck = () -> (
    if bmFamily == "sat" then (
	(M, f) := (bmLookupSat bmCaseNm)();
	a := msolveSaturate(bmCopy M, f);
	b := saturate(bmCopy M, f, Strategy => Iterate);
	g := msolveSaturateVerify(bmCopy M, f);
	bmSay(a == b, "guard=" | (if g === null then "declined"
		else if g == b then "accepted-correct" else "ACCEPTED-WRONG")))
    else if bmFamily == "elim" then (
	(bld, k) := toSequence bmElimCases#bmCaseNm;
	(R, I) := bld();
	ma := msolveEliminationGB(ideal generators I, take(gens R, k));
	mb := eliminate(take(gens R, k), ideal generators I);
	if ma === null then (bmSay(true, "msolve-declined"); return);
	-- msolveEliminationGB returns its answer in ring I already; the name map
	-- is only needed if some path ever returns it in the block order ring
	if ring ma =!= R then ma = (bmByName(R, ring ma)) ma;
	bmSay(ma == mb,
	    "msolve=" | toString numgens ma | " m2=" | toString numgens mb))
    else if bmFamily == "melim" then (
	(R3, m, k3) := (bmLookupMelim bmCaseNm)();
	ea := bmMelimMethods#"msolve" (R3, m, k3);
	eb := bmMelimMethods#"m2" (R3, m, k3);
	if ea === null then (bmSay(true, "msolve-declined"); return);
	eb = substitute(eb, vars R3);  -- same variables in the same order
	-- an element of the elimination part involves none of the eliminated
	-- variables, so both answers live over the subring; compare there
	elimv := set take(gens R3, k3);
	bad := select(flatten entries (ea | eb), p -> #(set support p * elimv) > 0);
	if #bad > 0 then (bmSay(false, "entry uses an eliminated variable"); return);
	B := (coefficientRing R3)(monoid [take(gens R3, {k3, numgens R3 - 1})]);
	psi := map(B, R3, join(apply(k3, i -> 0_B), apply(numgens B, i -> B_i)));
	bmSay(image psi ea == image psi eb,
	    "msolve=" | toString numcols ea | " m2=" | toString numcols eb))
    else error("unknown family " | bmFamily))

bmCheck()
exit 0
