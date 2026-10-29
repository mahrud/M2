-- Capture every matrix handed to syz during the pushforward, then time the
-- pieces of msolveSyzygyHook separately: msolve's own syzygy computation, and
-- the mingens minimalization the hook runs on top of msolve's output.
needsPackage "DirectSummands"
needsPackage "Msolve"
debug Msolve
p0 = value getenv "CHAR"
useM = (getenv "USE_MSOLVE") == "1"
if useM then setupMsolve()
capt = new MutableList from {}
oldSyz = lookup(syz, Matrix);
syz Matrix := opts -> m -> (capt#(#capt) = m; (oldSyz opts) m)
kk = ZZ/p0
R = quotient Grassmannian(1, 3, CoefficientRing => kk)
X = Proj R
t = elapsedTiming (M = module frobeniusPushforward(1, OO_X));
<< "##### char=" << p0 << " msolve=" << useM << " total=" << t#0
   << " nsyz=" << #capt << endl << flush
for i to #capt - 1 do (
    m := capt#i;
    if m == 0 or not isHomogeneous m or not msolveSyzygyUsable m
    then (<< "  mat " << i << " skipped" << endl; continue);
    tn := elapsedTiming (N := (oldSyz (options syz)) m);
    tm := elapsedTiming (S := msolveSyzygy m);
    tg := elapsedTiming (S2 := mingens image S);
    << "  mat " << i << " " << numrows m << "x" << numcols m
       << " | m2syz=" << tn#0 << " ncols=" << numcols N
       << " | msolveSyz=" << tm#0 << " ncols=" << numcols S
       << " | mingens=" << tg#0 << " ncols=" << numcols S2 << endl << flush)
exit 0
