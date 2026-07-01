debug Core
needsPackage "Polyhedra"

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
r = 2 -- number of parameters (= numcols ctx - 2)

-- rawBarvinokEnumerate returns a single ZZ matrix of tagged rows
-- [tag, chamberIndex, index, data...], entirely numeric (no strings):
--   tag 0 (facet):  data = facet coefficients (r values)
--   tag 1 (divdef): data = [coeff_0..coeff_{r-1}, const, denom] for the
--                   index-th floor()/div term's affine argument
--   tag 2 (term):   data = [coeffNum, coeffDen, exp_0..exp_{r-1},
--                   divExp_0..divExp_{maxDivs-1}]
ret = entries map(ZZ, rawBarvinokEnumerate(raw m, raw ctx))
maxDivs = (#ret#0) - 5 - r
assert(#unique apply(ret, row -> row#1) == 2)

-- reconstruct each chamber's cone and quasipolynomial function purely from
-- the numeric row data
chambers = hashTable apply(unique apply(ret, row -> row#1), ci -> (
	chRows := select(ret, row -> row#1 == ci);
	facetRows  := select(chRows, row -> row#0 == 0);
	divdefRows := sort select(chRows, row -> row#0 == 1);
	termRows   := select(chRows, row -> row#0 == 2);
	C := coneFromHData matrix apply(facetRows, row -> row_{3..3+r-1});
	divdefs := apply(divdefRows, row -> (toList row_{3..3+r-1}, row#(3+r), row#(3+r+1)));
	terms   := apply(termRows,   row -> (row#3, row#4, toList row_{5..5+r-1}, toList row_{5+r..5+r+maxDivs-1}));
	qpf := args -> (
	    ts := toList toSequence args;
	    sum(terms, (cNum, cDen, expL, divExpL) -> (
		    val := (cNum/cDen) * product(r, i -> ts#i^(expL#i));
		    val * product(#divExpL, j -> (
			    if divExpL#j == 0 then 1 else (
				(coefL, cst, denom) := divdefs#j;
				floor((sum(r, i -> coefL#i*ts#i) + cst)/denom)^(divExpL#j))))
		    )));
	ci => (C, qpf)
	))

-- brute force vector partition count, to check against the quasipolynomials
count = (t0, t1) -> sum(0..t1, v3 -> (
        v1max := t0 + 3*v3;
        if v1max < 0 then 0 else v1max + 1))

scan(values chambers, (C, f) ->
    scan({(4,1), (10,2), (7,3)}, (t0,t1) ->
        if contains(C, matrix{{t0},{t1}}) then
        assert(f(t0,t1) == count(t0,t1))))
