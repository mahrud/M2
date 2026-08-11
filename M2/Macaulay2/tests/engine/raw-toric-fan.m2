-- Engine-level fan refinement tests.

debug Core;

R = raw matrix {{1,0,1}, {0,1,1}, {-1,0,1}, {0,-1,1}};
C = {1,4,0,1,2,3};
rayRowsR = entries R;

S = map(ZZ, rawSimplicialFan(R, C, 0, 0, 0, 1, false));
assert(numRows S >= 6);
assert(all(0..3, i -> S_(i,0) == 0 and S_(i,1) == 3));
assert(all(0..3, i -> take(drop((entries S)#i, 2), 3) == rayRowsR#i));
assert(all(select(entries S, r -> r#0 == 1), r -> r#1 == 3));

T = map(ZZ, rawSmoothFan(R, C, 0, 0, 0, 1, false));
assert(numRows T >= numRows S);
assert(all(select(entries T, r -> r#0 == 1), r -> r#1 == 3));
rayRows = select(entries T, r -> r#0 == 0);
assert(all(0..3, i -> take(drop(rayRows#i, 2), 3) == rayRowsR#i));
assert(all(select(entries T, r -> r#0 == 1), r -> (
    coneRays := apply(take(drop(r, 2), r#1), i ->
        take(drop(rayRows#i, 2), (rayRows#i)#1));
    D := matrix coneRays;
    det D == 1 or det D == -1
    )));

U = map(ZZ, rawSmoothFan(R, C, 1, 17, 0, 1, false));
assert(all(select(entries U, r -> r#0 == 1), r -> r#1 == 3));
rayRowsU = select(entries U, r -> r#0 == 0);
assert(all(select(entries U, r -> r#0 == 1), r -> (
    D := matrix apply(take(drop(r, 2), r#1), i ->
        take(drop(rayRowsU#i, 2), (rayRowsU#i)#1));
    det D == 1 or det D == -1
    )));
