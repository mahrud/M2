-- http://nemocas.org/benchmarks.html
R = GF(17^11)
S = R[y]
T = S/(y^3 + 3*a*y + 1)
U = T[z]
f = (3*y^2 + y + a)*z^2 + ((a + 2)*y^2 + a + 1)*z + (4*a*y + 3)
g = (7*y^2 - y + 2*a + 7)*z^2 + (3*y^2 + 4*a + 1)*z + ((2*a + 1)*y + 1)
s = f^12
t = (s + g)^12
time r = resultant(s, t, y)
