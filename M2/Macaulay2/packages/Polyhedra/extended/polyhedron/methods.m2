 -- PURPOSE : Computing the smallest face of 'P' containing 'p'
--   INPUT : '(p,P)',  where 'p' is a point given as a matrix and
--     	    	       'P' is a polyhedron
--  OUTPUT : The smallest face containing 'p' as a polyhedron
smallestFace = method()
smallestFace(Matrix,Polyhedron) := (p,P) -> (
     -- Checking for input errors
     if numColumns p =!= 1 or numRows p =!= ambDim(P) then error("The point must lie in the same space");
     -- p = chkZZQQ(p,"point");
     -- Checking if 'P' contains 'p' at all
     if contains(P,convexHull p) then (
	      (M,v) := halfspaces P;
     	   (N,w) := hyperplanes P;
     	  -- Selecting the half-spaces that fullfil equality for p
	  -- and adding them to the hyperplanes
	  v = promote(v,QQ);
	  pos := select(toList(0..(numRows M)-1), i -> (M^{i})*p == v^{i});
	  N = N || M^pos;
	  w = w || lift(v^pos,ZZ);
	  polyhedronFromHData(M,lift(v,ZZ),N,w))
     else emptyPolyhedron ambDim(P))



-- PURPOSE : Checks if the polytope is normal
--   INPUT : 'P'  a Polyhedron, which must be compact
--  OUTPUT : 'true' or 'false'
-- COMMENT : The polytope is normal if the lattice of the cone over the polytope embedded on height 1 
--     	     is generated by the lattice points on height 1
isNormal Polyhedron := P -> getProperty(P, computedNormal)


-- PURPOSE : Computing the Ehrhart polynomial of a polytope
--   INPUT : 'P',  a polyhedron which must be compact, i.e. a polytope
--  OUTPUT : A polynomial in QQ[x], the Ehrhart polynomial
-- COMMENT : Compactness is checked within latticePoints
ehrhart = method(TypicalValue => RingElement)
ehrhart Polyhedron := P -> getProperty(P, computedEhrhart)


-- PURPOSE : Computing the polar of a given polyhedron
--   INPUT : 'P',  a Polyhedron
--  OUTPUT : A Polyhedron, the set { v | v*p<=1 forall p in P}
polar = method(TypicalValue => Polyhedron)
polar Polyhedron := P -> getProperty(P, computedPolar)


-- PURPOSE : Checks if a polytope is very ample
--   INPUT : 'P'  a Polyhedron, which must be compact
--  OUTPUT : 'true' or 'false'
isVeryAmple = method()
isVeryAmple Polyhedron := P -> getProperty(P, computedVeryAmple)


-- PURPOSE : Computing the vertex-edge-matrix of a polyhedron
--   INPUT : 'P',  a polyhedron
--  OUTPUT : a matrix, where the columns are indexed by the edges and the rows indexed by the vertices and has 1 as entry
--           if the corresponding edge contains this vertex
vertexEdgeMatrix = method(TypicalValue => Matrix)
vertexEdgeMatrix Polyhedron := P -> (
   -- list the edges and the vertices
   eP := apply(faces(dim P -1,P), f -> f#0);
   nEdge := #eP;
   nVert := numColumns vertices P;
   result := map (ZZ^nVert, ZZ^nEdge, 0);
   result = (matrix {{1..nEdge}}) || result;
   result = (transpose matrix {{0..nVert}}) | result;
   result = mutableMatrix result;
   i := 0;
   for edge in eP do (
      s := edge#0;
      t := edge#1;
      result_(s + 1, i + 1) = 1;
      result_(t + 1, i + 1) = 1;
      i = i+1;
   );
   matrix result
)


-- PURPOSE : Computing the vertex-facet-matrix of a polyhedron
--   INPUT : 'P',  a polyhedron
--  OUTPUT : a matrix, where the columns are indexed by the facets and the rows are indexed by the vertices and has 1 as entry
--           if the corresponding facet contains this vertex
vertexFacetMatrix = method(TypicalValue => Matrix)
vertexFacetMatrix Polyhedron := P -> (
   -- list the facets and the vertices
   fP := apply(faces(1,P), f -> f#0);
   nFacet := #fP;
   nVert := numColumns vertices P;
   result := map (ZZ^nVert, ZZ^nFacet, 0);
   result = (matrix {{1..nFacet}}) || result;
   result = (transpose matrix {{0..nVert}}) | result;
   result = mutableMatrix result;
   i := 0;
   for facet in fP do (
      for v in facet do (
         result_(v + 1, i + 1) = 1;
      );
      i = i+1;
   );
   matrix result
)


-- PURPOSE : Computing the face fan of a polytope
--   INPUT : 'P',  a Polyhedron, containing the origin in its interior
--  OUTPUT : The Fan generated by the cones over all facets of the polyhedron
faceFan = method(TypicalValue => Fan)
faceFan Polyhedron := P -> (
   -- Checking for input errors
   if not inInterior(map(QQ^(ambDim P),QQ^1,0),P) then  error("The origin must be an interior point.");
   if not isCompact P then error("Polyhedron must be compact");
   resultRays := vertices P;
   resultCones := faces(1, P);
   resultCones = apply(resultCones, c -> c#0);
   fan(resultRays, resultCones)
)


-- PURPOSE : Computing the cell decomposition of a compact polyhedron given by a weight vector on the lattice points
--   INPUT : '(P,w)',  where 'P' is a compact polyhedron and 'w' is a one row matrix with lattice points of 'P' 
--     	    	       many entries
--  OUTPUT : A list of polyhedra that are the corresponding cell decomposition
regularSubdivision = method(TypicalValue => List)
regularSubdivision (Polyhedron,Matrix) := (P,w) -> (
   if not isCompact P then error("The polyhedron must be compact.");
   n := dim P;
   LP := latticePoints P;
   LP = transpose matrix apply(LP, l -> flatten entries l);
   -- Checking for input errors
   if numColumns w != numColumns LP or numRows w != 1 then error("The weight must be a one row matrix with number of lattice points many entries");
   S := regularSubdivision (LP, w);
   apply (S, s -> convexHull LP_s)
)

regularSubdivision (Matrix,Matrix) := (MM, w) -> (
   M := promote(MM, QQ);
   n := numColumns M;
   -- Checking for input errors
   if numColumns w != numColumns M or numRows w != 1 then error("The weight must be a one row matrix with number of points many entries");
   P := convexHull(M||w,matrix (toList(numRows M:{0})|{{1}}));
   F := select(faces (1,P), f -> #(f#1) ==0);
   pointIndices := new MutableHashTable;
   for i from 0 to numcols M-1 do
      pointIndices#(M_{i}) = i;
   permutation := new MutableHashTable;
   vertP := (vertices P)^{0..numrows M - 1};
   for i from 0 to numcols vertP - 1 do
      permutation#i = pointIndices#(vertP_{i});
   sort apply (F, f -> sort apply(f#0, i->permutation#i))
  )




--   INPUT : 'P',  a polyhedron,
--  OUTPUT : A matrix, a basis of the sublattice spanned by the lattice points of 'P'
sublatticeBasis Polyhedron := P -> (
     L := latticePoints P;
     -- Checking for input errors
     if L == {} then error("The polytope must contain lattice points.");
     -- Translating 'P' so that it contains the origin if it did not already
     if all(L,l -> l != 0) then L = apply(L, l -> l - L#0);
     sublatticeBasis(matrix {L}))
   
-- PURPOSE : Calculating the preimage of a polytope in the sublattice generated by its lattice points
--   INPUT : 'P',  a polyhedron
--  OUTPUT : A polyhedron, the projected polyhedron, which is now normal
toSublattice = method()
toSublattice Polyhedron := P -> (
     L := latticePoints P;
     -- Checking for input errors
     if L == {} then error("The polytope must contain lattice points.");
     b := L#0;
     -- Translating 'P' so that it contains the origin if it did not already
     if all(L,l -> l != 0) then L = apply(L, l -> l - L#0);     
     affinePreimage(sublatticeBasis matrix {L},P,b))


-- PURPOSE : Compute the corresponding face of the polar polytope
--   INPUT : 'P',  a Polyhedron
--  OUTPUT : A Polyhedron, if 'P' is the face of some polyhedron 'Q' then the
--     	     result is the dual face on the polar of 'Q'. If 'P' is not a face
--           then it is considered as the face of itself and thus the 
--           polarFace is the empty Polyhedron
polarFace = method(TypicalValue => Polyhedron)
polarFace(Polyhedron, Polyhedron) := (f, P) -> (
   if not isFace(f, P) then error("First polyhedron should be face of second.");
   if f == P then return emptyPolyhedron dim f;
   --local faceOf;
   V := transpose vertices f;
   R := transpose rays f;
   Pd := polar P;
   codimensionfd := dim f - (numColumns linealitySpace P) + 1;
   L := facesAsPolyhedra(codimensionfd, Pd);
   Pd = first select(1,L, l -> all(flatten entries(V*(vertices l)),e -> e == -1) and V*(rays l) == 0 and R*(vertices l | rays l) == 0);
   Pd
)	       


-- PURPOSE : Checks if a lattice polytope is reflexive
--   INPUT : 'P'  a Polyhedron
--  OUTPUT : 'true' or 'false'
isReflexive = method(TypicalValue => Boolean)
isReflexive Polyhedron := (cacheValue symbol isReflexive)(P -> isLatticePolytope P and inInterior(matrix toList(ambDim P:{0}),P) and isLatticePolytope polar P)


-- PURPOSE : Triangulating a compact Polyhedron
--   INPUT : 'P',  a Polyhedron
--  OUTPUT : A list of the simplices of the triangulation. Each simplex is given by a list 
--    	     of its vertices.
--COMMENTS : The triangulation is build recursively, for each face that is not a simplex it takes 
--     	     the weighted centre of the face. for each codim 1 face of this face it either takes the 
--     	     convex hull with the centre if it is a simplex or triangulates this in the same way.
barycentricTriangulation = method()
barycentricTriangulation Polyhedron := P -> (
     -- Defining the recursive face triangulation
     -- This takes a polytope and computes all facets. For each facet that is not a simplex, it calls itself
     -- again to replace this facet by a triangulation of it. then it has a list of simplices triangulating 
     -- the facets. Then it computes for each of these simplices the convex hull with the weighted centre of 
     -- the input polytope. The weighted centre is the sum of the vertices divided by the number of vertices.
     -- It returns the resulting list of simplices in a list, where each simplex is given by a list of its 
     -- vertices.
     -- The function also needs the dimension of the Polyhedron 'd', the list of facets of the original 
     -- polytope, the list 'L' of triangulations computed so far and the dimension of the original Polytope.
     -- 'L' contains a hash table for each dimension of faces of the original Polytope (i.e. from 0 to 'n').
     -- If a face has been triangulated than the list of simplices is saved in the hash table of the 
     -- corresponding dimension with the weighted centre of the original face as key.
     recursiveFaceTriangulation := (P,d,originalFacets,L,n) -> (
	  -- Computing the facets of P, given as lists of their vertices
	  F := intersectionWithFacets({(set P,set{})},originalFacets);
	  F = apply(F, f -> toList(f#0));
	  d = d-1;
	  -- if the facets are at least 2 dimensional, then check if they are simplices, if not call this 
	  -- function again
	  if d > 1 then (
	       F = flatten apply(F, f -> (
			 -- Check if the face is a simplex
			 if #f != d+1 then (
			      -- Computing the weighted centre
			      p := (sum f)*(1/#f);
			      -- Taking the hash table of the corresponding dimension
			      -- Checking if the triangulation has been computed already
			      if L#d#?p then L#d#p
			      else (
				   -- if not, call this function again for 'f' and then save this in 'L'
				   (f,L) = recursiveFaceTriangulation(f,d,originalFacets,L,n);
				   L = merge(L,hashTable {d => hashTable{p => f}},(x,y) -> merge(x,y,));
				   f))
			 else {f})));
	  -- Adding the weighted centre to each face simplex
	  q := (sum P)*(1/#P);
	  P = apply(F, f -> f | {q});
	  (P,L));
     -- Checking for input errors
     if not isCompact P then error("The polytope must be compact!");
     n := dim P;
     -- Computing the facets of P as lists of their vertices
     (HS,v) := halfspaces P;
     (hyperplanesTmp,w) := hyperplanes P;
     originalFacets := apply(numRows HS, i -> polyhedronFromHData(HS,v, hyperplanesTmp || HS^{i}, w || v^{i}));
     originalFacets = apply(originalFacets, f -> (
	       V := vertices f;
	       (set apply(numColumns V, i -> V_{i}),set {})));
     -- Making a list of the vertices of P
     P = vertices P;
     P = apply(numColumns P, i -> P_{i});
     if #P == n+1 then {P} else (
	  d := n;
	  -- Initiating the list of already computed triangulations
	  L := hashTable apply(n+1, i -> i => hashTable {});
	  (P,L) = recursiveFaceTriangulation(P,d,originalFacets,L,n);
	  P))


latticeVolume Polyhedron := P -> getProperty(P, latticeVolume)


-- PURPOSE : Computing the volume of a full dimensional polytope
--   INPUT : 'P',  a compact polyhedron
--  OUTPUT : QQ, giving the volume of the polytope
volume = method(TypicalValue => QQ)
volume Polyhedron := P -> (
   if isEmpty P then return 0;
   d := dim P;
   result := latticeVolume P;
   result/(d!)
)


	       
-- PURPOSE : Computing a point in the relative interior of a Polyhedron 
--   INPUT : 'P',  a Polyhedron
--  OUTPUT : 'p',  a point given as a matrix
interiorPoint = method(TypicalValue => Matrix)
interiorPoint Polyhedron := P -> (
     -- Checking for input errors
     if isEmpty P then error("The polyhedron must not be empty");
     Vm := vertices P | promote(rays P,QQ);
     n := numColumns Vm;
     ones := matrix toList(n:{1/n});
     -- Take the '1/n' weighted sum of the vertices
     Vm * ones)


-- PURPOSE : Computing the face of a Polyhedron where a given weight attains its minimum
--   INPUT : '(v,P)',  a weight vector 'v' given by a one column matrix over ZZ or QQ and a 
--     	     	       Polyhedron 'P'
--  OUTPUT : a Polyhedron, the face of 'P' where 'v' attains its minimum
minFace (Matrix,Polyhedron) := (v,P) -> (
   -- Checking for input errors
   if numColumns v =!= 1 or numRows v =!= ambDim(P) then error("The vector must lie in the same space as the polyhedron");
   C := dualCone tailCone P;
   V := vertices P;
   R := rays P;
   LS := linSpace P;
   -- The weight must lie in the dual of the tailcone of the polyhedron, otherwise there is 
   -- no minimum and the result is the empty polyhedron
   if contains(C,v) then (
      -- Compute the values of 'v' on the vertices of 'V'
      Vind := flatten entries ((transpose v)*V);
      -- Take the minimal value(s)
      Vmin := min Vind;
      Vind = positions(Vind, e -> e == Vmin);
      -- If 'v' is in the interior of the dual tailCone then the face is exactly spanned 
      -- by these vertices
      if inInterior(v,C) then convexHull(V_Vind,LS | -LS)
      else (
         -- Otherwise, one has to add the rays of the tail cone that are orthogonal to 'v'
         Rind := flatten entries ((transpose v)*R);
         Rind = positions(Rind, e -> e == 0);
         convexHull(V_Vind,R_Rind | LS | -LS)
      )
   ) else emptyPolyhedron ambDim P
)


-- PURPOSE : Computing the face of a Polyhedron where a given weight attains its maximum
--   INPUT : '(v,P)',  a weight vector 'v' given by a one column matrix over ZZ or QQ and a 
--     	     	       Polyhedron 'P'
--  OUTPUT : a Polyhedron, the face of 'P' where 'v' attains its maximum
maxFace (Matrix,Polyhedron) := (v,P) -> minFace(-v,P)


--   INPUT : '(p,P)',  where 'p' is a point given by a matrix and 'P' is a Polyhedron
--  OUTPUT : 'true' or 'false'
inInterior (Matrix,Polyhedron) := (p,P) -> (
     hyperplanesTmp := hyperplanes P;
     hyperplanesTmp = (hyperplanesTmp#0 * p)-hyperplanesTmp#1;
     all(flatten entries hyperplanesTmp, e -> e == 0) and (
	  HS := halfspaces P;
	  HS = (HS#0 * p)-HS#1;
	  all(flatten entries HS, e -> e < 0)))



-- PURPOSE : Computing the affine hull
--   INPUT : 'P',  a Polyhedron
--  OUTPUT : the affine hull of 'P' as a Polyhedron
affineHull = method(TypicalValue => Polyhedron)
affineHull Polyhedron := P -> (
     M := vertices P;
     R := promote(rays P,QQ);
     -- subtracting the first vertex from all other vertices
     N := (M+M_{0}*(matrix {toList(numColumns M:-1)}))_{1..(numColumns M)-1};
     convexHull(M_{0},N | -N | R | -R));
 

