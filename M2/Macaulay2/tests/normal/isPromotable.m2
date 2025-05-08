debug Core
assert all({ZZ, ZZ/11, QQ, RR_53, CC_53, GF 4},
    k -> isPromotable(k, k[x,y]) and not isPromotable(k[x,y], k))

-- from EngineTests/Ring.Test.RR.CC.m2
RCC = CC_53[x,y,z];
p = 0_(coefficientRing RCC)
assert(promote(p, RCC) == 0_RCC)
