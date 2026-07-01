debug Core

-- basis polyhedron for a vector partition function of 4 generators mapped to
-- 2 parameters via v0+v1-3v3=t0, v2+v3=t1 (BarvinokFeature.md's worked example)
m = matrix {
  {0, 1, 1, 0, -3, -1, 0, 0},
  {0, 0, 0, 1,  1,  0,-1, 0},
  {1, 1, 0, 0,  0,  0, 0, 0},
  {1, 0, 1, 0,  0,  0, 0, 0},
  {1, 0, 0, 1,  0,  0, 0, 0},
  {1, 0, 0, 0,  1,  0, 0, 0}
  }
ctx = map(ZZ^0, ZZ^4, 0)

ret = rawBarvinokEnumerate(raw m, raw ctx)
assert(#ret == 2)

-- brute force vector partition count, to check against the quasipolynomials
count = (t0, t1) -> sum(0..t1, v3 -> (
        v1max := t0 + 3*v3;
        if v1max < 0 then 0 else v1max + 1))

needsPackage "Polyhedra"
scan(ret, (facets, str) -> (
    C := coneFromHData map(ZZ, facets);
    f := (value str) @@ toSequence;
    scan({(4,1), (10,2), (7,3)}, (t0,t1) ->
        if contains(C, matrix{{t0},{t1}}) then
        assert(f(t0,t1) == count(t0,t1)))
    ))
