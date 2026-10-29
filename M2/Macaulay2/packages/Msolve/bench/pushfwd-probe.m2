needsPackage "DirectSummands"
needsPackage "Msolve"
debugLevel = 2
p = value getenv "CHAR"
useM = (getenv "USE_MSOLVE") == "1"
if useM then setupMsolve()
<< "##### char=" << p << " msolve=" << useM << endl << flush
kk = ZZ/p
R = quotient Grassmannian(1, 3, CoefficientRing => kk)
X = Proj R
t = elapsedTiming (M = module frobeniusPushforward(1, OO_X));
<< "##### TOTAL char=" << p << " msolve=" << useM << " secs=" << t#0 << endl << flush
<< "##### rank=" << rank M << endl << flush
exit 0
