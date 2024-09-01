needs "methods.m2"

-----------------------------------------------------------------------------
-- AtomicInt
-----------------------------------------------------------------------------

AtomicInt.synonym = "atomic integer"

scan({symbol +=, symbol -=, symbol &=, symbol |=, symbol ^^=},
    op -> typicalValues#(op, AtomicInt) = ZZ)

store = method()
store(AtomicInt, ZZ) := atomicStore

exchange = method()
exchange(AtomicInt, ZZ) := atomicExchange

compareExchange = method()
compareExchange(AtomicInt, ZZ, ZZ) := atomicCompareExchange

-----------------------------------------------------------------------------
-- Mutex
-----------------------------------------------------------------------------

Mutex.synonym = "mutex"
globalAssignment Mutex
net Mutex := x -> toString (
    if hasAttribute(x, ReverseDictionary)
    then getAttribute(x, ReverseDictionary)
    else x)

lock = method()
lock Mutex := M -> ( lock0 M; M )
lock(Mutex, Function) := lockFunction -- defined in pthread.d
lock Function := f -> lockFunction(new Mutex, f)

tryLock = method()
tryLock Mutex := M -> ( tryLock0 M; M )

unlock = method()
unlock Mutex := unlock0

-----------------------------------------------------------------------------
-- async/await
-----------------------------------------------------------------------------
-- TODO: unreachable tasks should be cancelled

-- TODO: move these to the interpreter
async = method(Dispatch => Thing)
async Function := Function => f -> x -> schedule(f, x)

-- TODO: what other data structures can hold tasks?
await = method(Dispatch => Thing)
await Task      := Thing     => await @@ taskResult
await Type      :=
await Thing     := Thing     => identity
await BasicList := BasicList => L -> apply(L, await)
-- TODO: should this also replace the keys and handle conflicts?
await HashTable := HashTable => H -> (
    if class H =!= HashTable then H else applyValues(H, await))
-- TODO: applyValues should take a mutable hash table, and perhaps modify it in place?
await MutableHashTable := MutableHashTable => H -> (
    if class H =!= MutableHashTable then H else ( scan(keys H, k -> H#k = await H#k); H ))

-----------------------------------------------------------------------------

parallelApplyRaw = (L, f) ->
     -- 'reverse's to minimize thread switching in 'taskResult's:
     reverse (taskResult \ reverse apply(L, e -> schedule(f, e)));
parallelApply = method(Options => {Strategy => null})
parallelApply(BasicList, Function) := o -> (L, f) -> (
     if o.Strategy === "raw" then return parallelApplyRaw(L, f);
     n := #L;
     numThreads := min(n + 1, maxAllowableThreads);
     oldAllowableThreads := allowableThreads;
     if allowableThreads < numThreads then allowableThreads = numThreads;
     numChunks := 3 * numThreads;
     res := if n <= numChunks then toList parallelApplyRaw(L, f) else
	  flatten parallelApplyRaw(pack(L, ceiling(n / numChunks)), chunk -> apply(chunk, f));
     allowableThreads = oldAllowableThreads;
     res);
