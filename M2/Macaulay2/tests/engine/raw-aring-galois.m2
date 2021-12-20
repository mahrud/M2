debug Core
p = 3
d = 2
R = ZZ/p[x]
A = rawARingGaloisField(p, d)
c = rawARingGFPolynomial A
f = sum(#c, i -> c#i * x^i)
a = rawRingVar(A, 0)

rawARingGFCoefficients(3*a^3 + 2*a^2)
4_A
A_0

--------
A = quotient ideal conwayPolynomial(3,2)
R = rawARingGaloisFieldM2 raw A_0
a = rawRingVar(R, 0)

--------
R1 = ZZ/2[x]/ideal"x8+x4+x3+x1+1"
assert     isPrimitive(x+1)
assert not isPrimitive x
--findPrimitive R1

R2 = ZZ/2[x]/ideal"x8+x4+x3+x2+1"
assert not isPrimitive(x+1)
assert     isPrimitive x
--findPrimitive R1

rawR = rawARingGaloisFieldM2 raw(x_R2)
rawR = rawARingGaloisFieldFlintBig raw(x_R2)
rawR = rawARingGaloisFieldFlintZech raw(x_R2)
rawR = rawARingGaloisField(2, 8)

findPrimitive baseRing GF(2, 8)

end--
-- get these two to work
A1 = GF R1
assert     isPrimitive(x+1)
assert not isPrimitive x
findPrimitive R1

A2 = GF R2
assert not isPrimitive(x+1)
assert     isPrimitive x
findPrimitive R1
