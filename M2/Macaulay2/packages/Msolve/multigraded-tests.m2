-- Checks for the multigraded Betti and Hilbert entry points.
--
-- These exercise the binding rather than the engine: msolve's own selftest
-- pins down the tables themselves, and what is checked here is that the ring's
-- grading is marshalled across intact and that what comes back lands on
-- Macaulay2's scale.  Kept as a script rather than folded into tests.m2 so the
-- slow toric example at the end can be skipped.
--
-- Run with: M2 --script multigraded-tests.m2

debug Core
needsPackage "Msolve"
errorDepth = 2

nfail = 0
check1 = (name, val) -> (
    if val then print("ok   " | name)
    else (nfail = nfail + 1; print("FAIL " | name)))

-- There is one table and one entry point: rawMsolveMinimalBetti reports the
-- multidegrees themselves, a singly graded ring being the case where those are
-- vectors of length one, so a key is always (level, multidegree, heft degree).
-- What comes back is a plain BettiTally either way; printing such a table with
-- its multidegrees rather than folded onto the heft degree is `multigraded`,
-- which is the caller's to ask for.
mgBetti = (m, len) -> unpackMsolveBetti rawMsolveMinimalBetti(raw m, len, 1, 0)

-- summing that table over each heft fibre is what the heft indexed table is,
-- the one Macaulay2's own minimalBetti reports for a multigraded module
collapse = B -> new BettiTally from applyKeys(B, k -> (k#0, {k#2}, k#2), plus)
hBetti = (m, len) -> collapse mgBetti(m, len)

-- msolve puts the heft degree in the {d} slot of a collapsed key, which is
-- where Macaulay2's own minimalBetti puts it too; M2 fills the third slot with
-- a Weights scaled number instead, so the two are compared on (level, degree).
strip = B -> new BettiTally from apply(pairs B, (k, v) -> (k#0, k#1) => v)

-- the multigraded numerator is the alternating sum of the multigraded table
num = (B, R) -> sum(pairs B, (k, v) -> (-1)^(k#0) * v * product(#k#1, i -> R_i^(k#1#i)))

-----------------------------------------------------------------------------
printerr "-- 1. standard graded, the pre-existing case"
S = ZZ/32003[x,y,z,w]
A = matrix {{x,y,z}, {y,z,w}}
check1("standard betti",    mgBetti(A, 5) === minimalBetti coker A)
check1("standard poincare", rawMsolvePoincare(raw A, 1, 0) === raw poincare coker A)
-- a singly graded ring is the r = 1 case: the degree is the vector {d}, which
-- is what the collapsed key holds too, so collapsing is the identity here
check1("standard collapse", hBetti(A, 5) === mgBetti(A, 5))
check1("standard is plain", not instance(mgBetti(A, 5), MultigradedBettiTally))

-----------------------------------------------------------------------------
printerr "-- 2. singly graded but weighted: degrees agree with the grevlex weights"
Rw = ZZ/32003[a,b,c, Degrees => {1,2,3}]
Iw = ideal(a*c - b^2, b^3 - c^2)
check1("weighted betti",    hBetti(gens Iw, 5) === minimalBetti coker gens Iw)
check1("weighted poincare", rawMsolvePoincare(raw gens Iw, 1, 0) === raw poincare coker gens Iw)

-----------------------------------------------------------------------------
printerr "-- 3. weights unrelated to the grading: refused before, allowed now"
Ru = ZZ/32003[a,b,c, Degrees => {1,2,3}, MonomialOrder => {GRevLex => {1,1,1}}]
Iu = ideal(a*c - b^2, b^3 - c^2)
check1("unrelated betti",    hBetti(gens Iu, 5) === minimalBetti coker gens Iu)
check1("unrelated poincare", rawMsolvePoincare(raw gens Iu, 1, 0) === raw poincare coker gens Iu)

-----------------------------------------------------------------------------
printerr "-- 4. P^1 x P^1, the smallest multigraded ring"
T = ZZ/32003[a,b,c,d, Degrees => {{1,0},{1,0},{0,1},{0,1}}]
J = ideal(a*c, b*d, a*d)
check1("P1xP1 betti",    strip hBetti(gens J, 5) === strip minimalBetti coker gens J)
check1("P1xP1 poincare", rawMsolvePoincare(raw gens J, 1, 0) === raw poincare coker gens J)
MB = mgBetti(gens J, 5)
-- the table is plain whatever the grading, and asking for the net that keeps
-- the multidegrees rather than folding them onto the heft degree is the
-- caller's `multigraded`, as it is for any other BettiTally
check1("P1xP1 is plain",       not instance(MB, MultigradedBettiTally))
check1("P1xP1 is multigraded", instance(multigraded MB, MultigradedBettiTally))
check1("P1xP1 numerator", num(MB, degreesRing T) == poincare coker gens J)
-- The load bearing check that the table really is multigraded rather than a
-- singly graded one wearing multidegrees: an actual minimal free resolution,
-- compared entry by entry.  Note {1,2} and {2,1} both sit over heft degree 3,
-- so a table that had collapsed onto the heft degree could not reproduce this.
check1("P1xP1 vs resolution", MB === betti freeResolution coker gens J)

-- a multigraded module with shifted, non-equal row degrees
N = coker matrix{{a*c, b*d}}
check1("module betti",     strip hBetti(presentation N, 5) === strip minimalBetti N)
check1("module poincare",  rawMsolvePoincare(raw presentation N, 1, 0) === raw poincare N)
MB = mgBetti(presentation N, 5)
check1("module numerator", num(MB, degreesRing T) == poincare N)
check1("module vs resolution", MB === betti freeResolution N)

-----------------------------------------------------------------------------
printerr "-- 5. a multigraded free module: the zero submodule short circuit"
F = T^{{-1,0},{0,-2},{-1,0}}
Z = map(F, T^0, 0)
check1("free betti",       strip hBetti(Z, 5) === strip minimalBetti image id_F)
check1("free poincare",    rawMsolvePoincare(raw Z, 1, 0) === raw poincare image id_F)
check1("free multigraded", mgBetti(Z, 5) === betti F)

-----------------------------------------------------------------------------
printerr "-- 6. a grading with no positive heft is refused, not answered"
G = ZZ/32003[u,v, Degrees => {1,-1}]
check1("no heft refused", null === try hBetti(vars G, 5) else null)

-----------------------------------------------------------------------------
printerr "-- 7. a three dimensional grading with a nontrivial heft"
X = ZZ/32003[p_0..p_4, Degrees => {{1,0,0},{1,0,0},{-2,1,0},{1,-1,1},{0,0,1}}]
K = ideal(p_0*p_2, p_1*p_3*p_4, p_0*p_1)
check1("heft3 betti",     strip hBetti(gens K, 6) === strip minimalBetti coker gens K)
check1("heft3 poincare",  rawMsolvePoincare(raw gens K, 1, 0) === raw poincare coker gens K)
check1("heft3 vs resolution", mgBetti(gens K, 0) === betti freeResolution coker gens K)

-----------------------------------------------------------------------------
printerr "-- 8. a truncated table is still exact as far as it goes"
check1("truncated",
    (new BettiTally from select(pairs mgBetti(gens J, 1), (k,v) -> k#0 <= 1))
    === new BettiTally from select(pairs mgBetti(gens J, 5), (k,v) -> k#0 <= 1))

-----------------------------------------------------------------------------
printerr "-- 9. the toric example from tests.m2 (the slow one)"
needsPackage "NormalToricVarieties"
Y = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/101)
SY = ring Y
M = truncate({4,4,4}, SY^2);
elapsedTime MB = mgBetti(presentation M, 8);
print multigraded MB
print peek MB
check1("toric betti",     strip collapse MB === strip minimalBetti M)
check1("toric poincare",  rawMsolvePoincare(raw presentation M, 3, 0) === raw poincare M)

print("")
print if nfail == 0 then "ALL CHECKS PASSED" else (toString nfail | " CHECKS FAILED")
