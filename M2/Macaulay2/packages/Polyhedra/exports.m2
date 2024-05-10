export {

   "PolyhedralObject", 
-- Cone object with associated methods:
   "Cone", 
   "coneFromVData",
   "coneFromRays" => "coneFromVData", 
   "coneFromHData",
   "coneFromInequalities" => "coneFromHData",
 
-- Polyhedron object with associated methods:
   "Polyhedron", 
   "convexHull",
   "polyhedron",
   "polyhedronFromHData",
   "polyhedronFromInequalities" => "polyhedronFromHData",
   "dualFaceRepresentationMap",
   "facets",
   "Fan", 
   "PolyhedralComplex", 
   "fan", 
   "fanFromGfan",
   "addCone", 
   "polyhedralComplex", 
   "addPolyhedron", 
   "ambDim", 
   "cones", 
   "maxCones", 
   "maxPolyhedra", 
   "halfspaces", 
   "hyperplanes", 
   "linSpace", 
   "linealitySpace",
   "linearTransform",
   "polyhedra", 
   "vertices",
   "areCompatible", 
   "commonFace", 
   "contains", 
   "isCompact", 
   "isComplete", 
   "isFace", 
   "isFullDimensional",
   "isLatticePolytope", 
   "isPointed", 
   "isPolytopal", 
   "isPure",
   "isReflexive",
   "isSimplicial", 
   --"isSmooth",
   "isVeryAmple",
   "faces", 
   "facesAsPolyhedra",
   "facesAsCones",
   "fVector", 
   "incompCones",
   "incompPolyhedra",
   "inInterior", 
   "interiorPoint", 
   "interiorVector",
   "interiorLatticePoints", 
   "latticePoints", 
   "latticeVolume",
   "maxFace", 
   "minFace", 
   "objectiveVector",
   "minkSummandCone", 
   "mixedVolume", 
   "polytope", 
   "posHull" => "coneFromVData",
   "proximum", 
   "skeleton", 
   "smallestFace", 
   "smoothSubfan", 
   "stellarSubdivision", 
   "tailCone", 
   "triangulate", 
   "volume", 
   "vertexEdgeMatrix", 
   "vertexFacetMatrix", 
   "affineHull", 
   "affineImage", 
   "affinePreimage", 
   "bipyramid", 
   "ccRefinement", 
   "directProduct",
   "dualCone", 
   "faceFan", 
   "imageFan", 
   "minkowskiSum", 
   "normalFan", 
   "nVertices",
   "polar",
   "polarFace", 
   "pyramid", 
   "sublatticeBasis", 
   "toSublattice",
   "crossPolytope", 
   "cellDecompose", 
   "cyclicPolytope", 
   "ehrhart", 
   "emptyPolyhedron", 
   "hirzebruch", 
   "hypercube", 
   "newtonPolytope", 
   "posOrthant", 
   "secondaryPolytope", 
   "statePolytope", 
   "stdSimplex",
   "simplex",
   "saveSession",
   "regularTriangulation",
   "barycentricTriangulation",
   "regularSubdivision",
   "minimalNonFaces",
   "stanleyReisnerRing"
}

