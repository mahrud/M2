--------------------------------
--- PDState commands -----------
--------------------------------
PDState = new Type of MutableHashTable
createPDState = method()
createPDState Ideal := I -> (
    new PDState from {
	"OriginalIdeal"     => I,
	"PrimesSoFar"       => new MutableHashTable from {},
	"IntersectionSoFar" => ideal (1_(ring I)),
	"isPrimeIdeal"      => true,
	"isPrimaryIdeal"    => true,
	"PrunedViaCodim"    => 0,
	-- mean wall time of the previous round of each strategy, used by
	-- splitIdeals to decide whether spawning tasks is worth the overhead
	"RoundCost"         => new MutableHashTable from {}
	}
    )

-- 'ideal p' (a chain of saturations) and 'codim p' (a Groebner basis) are the
-- expensive part of recording a prime, but they depend only on p and mutate only
-- p's own cache.  They are computed here, outside pdState, so that the worklist
-- loop in splitIdeals can do them in the per-ideal task and leave updatePDState
-- as a cheap, purely bookkeeping step run on the main thread.
preparePrimes = method()
preparePrimes List := L -> apply(L, p -> (p, ideal p, codim p))

updatePDState = method()
updatePDState (PDState,List,ZZ) := (pdState,L,pruned) -> (
  -- this function updates the pdState with the new primes in the list L, which
  -- consists of triples (p, ideal p, codim p) as returned by preparePrimes,
  -- all of which are known to be prime.
  -- Note: keep this free of engine calls; see preparePrimes above.
  ansSoFar := pdState#"PrimesSoFar";
  pdState#"PrunedViaCodim" = pdState#"PrunedViaCodim" + pruned;
  for t in L do (
     (p, I, c) := t;
     --pdState#"IntersectionSoFar" = trim intersect(pdState#"IntersectionSoFar", I);
     --<< endl << "   Adding codimension " << c << " prime ideal." << endl;
     if not ansSoFar#?c then
        ansSoFar#c = {(p,I)}
     else
        ansSoFar#c = append(ansSoFar#c,(p,I));
  );
  -*
  -- here we update the isPrime flag if L comes in with more than one
  -- prime, then the ideal is neither prime nor primary.
  -- the reason for this is that no single step will produce multiple redundant
  -- primes.  The only possible redundancy occurs when a prime is
  -- already in pdState and also comes into the list L
  if #L > 1 then (
     pdState#"isPrimeIdeal" = false;
     pdState#"isPrimaryIdeal" = false;
  );
  *-
)

numPrimesInPDState = method()
numPrimesInPDState PDState := pdState -> sum apply(pairs (pdState#"PrimesSoFar"), p -> #(p#1))

getPrimesInPDState = method()
getPrimesInPDState PDState := pdState ->
   flatten apply(pairs (pdState#"PrimesSoFar"), p -> (p#1) / last)

isRedundantIdeal = method()
isRedundantIdeal (AnnotatedIdeal,PDState) := (I,pdState) -> (
   -- Note: this is a read of pdState, and as written it is a no-op, so splitIdeals
   -- applies it before the (possibly parallel) round rather than inside it.  If the
   -- commented-out containment test below is ever restored, this becomes a read
   -- dependency on state that other branches are concurrently growing, and the
   -- worklist loop has to be revisited.
   -- the reason for this line is that once IndependentSet has been called, then
   -- I.Ideal no longer reflects the ideal
   if I.?IndependentSet then return false;
   primeList := getPrimesInPDState(pdState);
   -- as of now, all ideals are declared not redundant
   false
   -- this commented line takes too long!
   --any(primeList,p -> isSubset(p,I))
   )

flagPrimality = method()
flagPrimality(PDState, Boolean) := (pdState, primality) -> (pdState#"isPrime" = primality;)
