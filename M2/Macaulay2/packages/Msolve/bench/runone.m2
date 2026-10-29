-- One measurement per process: CASE and METHOD come from the environment, the
-- process prints a single CSV line and exits.  Running each measurement in its
-- own process is what keeps a Groebner basis cached by one bmMethod out of the
-- next bmMethod's timing, and lets the driver impose a hard timeout that an
-- in-process alarm cannot (msolve does not poll for interrupts inside its own
-- C loops).

needsPackage "Msolve"
debug Msolve  -- msolveSaturateVerify and msolveElimModuleGB are not exported

load "cases.m2"

bmFamily  = getenv "FAMILY"
bmCaseNm  = getenv "CASE"
bmMethod  = getenv "METHOD"
bmMode    = if getenv "MODE" == "" then "time" else getenv "MODE"

---------------------------------------------------------------------------

bmFingerprint = X -> (
    g := if instance(X, Module) then sort flatten entries generators gb X
    else if instance(X, Ideal) then sort flatten entries generators gb X
    else sort flatten entries X;
    a := ascii toString g;
    toString(#a) | ":" | toString sum(#a, i -> (i + 1) * a#i % 1073741827))

bmEmit = (status, secs, cpu, detail) -> (
    << "RESULT," << bmFamily << "," << bmCaseNm << "," << bmMethod << ","
    << bmNthread << "," << secs << "," << cpu << "," << status << "," << detail << endl;
    << flush)

---------------------------------------------------------------------------
bmRun = () -> (
    local val;
    if bmFamily == "sat" then (
	(M, f) := (bmLookupSat bmCaseNm)();
	strat := bmSatMethods#bmMethod;
	c0 := cpuTime();
	el := elapsedTiming (val = strat(M, f));
	c1 := cpuTime();
	if val === null then bmEmit("rejected", el#0, c1 - c0, "guard-declined")
	else bmEmit("ok", el#0, c1 - c0,
	    if bmMode == "check" then bmFingerprint val else toString numcols generators val))
    else if bmFamily == "elim" then (
	(bld, k) := toSequence bmElimCases#bmCaseNm;
	(R, I) := bld();
	strat2 := bmElimMethods#bmMethod;
	d0 := cpuTime();
	el2 := elapsedTiming (val = strat2(R, I, k));
	d1 := cpuTime();
	bmEmit("ok", el2#0, d1 - d0,
	    if bmMode == "check" then bmFingerprint val else toString numgens val))
    else if bmFamily == "melim" then (
	(R2, m, k2) := (bmLookupMelim bmCaseNm)();
	strat3 := bmMelimMethods#bmMethod;
	e0 := cpuTime();
	el3 := elapsedTiming (val = strat3(R2, m, k2));
	e1 := cpuTime();
	if val === null then bmEmit("rejected", el3#0, e1 - e0, "msolve-declined")
	else bmEmit("ok", el3#0, e1 - e0,
	    if bmMode == "check" then bmFingerprint val else toString numcols val))
    else error("unknown bmFamily " | bmFamily))

r = try bmRun() else bmEmit("error", -1, -1, "exception")
exit 0
