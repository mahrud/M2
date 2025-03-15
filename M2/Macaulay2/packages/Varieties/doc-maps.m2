-- (sheaf, Variety, Matrix)
-- (sheaf, Variety, Matrix, ZZ)
-- warning: missing node: sheaf(Matrix) cited by Matrix
-- warning: missing node: sheaf(Matrix,ZZ) cited by Matrix

-- warning: missing node: - SheafMap cited by -
-- warning: missing node: CoherentSheaf ** SheafMap cited by **
-- warning: missing node: Ext^ZZ(CoherentSheaf,SheafMap) cited by Matrix
-- warning: missing node: HH^ZZ SheafMap cited by Matrix
-- warning: missing node: RingElement * SheafMap cited by *
-- warning: missing node: SheafMap * RingElement cited by *
-- warning: missing node: SheafMap * SheafMap cited by *
-- warning: missing node: SheafMap * ZZ cited by *
-- warning: missing node: SheafMap ** CoherentSheaf cited by **
-- warning: missing node: SheafMap ** SheafMap cited by **
-- warning: missing node: SheafMap ** SheafOfRings cited by **
-- warning: missing node: SheafMap + SheafMap cited by +
-- warning: missing node: SheafMap ++ SheafMap cited by ++
-- warning: missing node: SheafMap - SheafMap cited by -
-- warning: missing node: SheafMap // SheafMap cited by //
-- warning: missing node: SheafMap == SheafMap cited by Boolean
-- warning: missing node: SheafMap == ZZ cited by Boolean
-- warning: missing node: SheafMap ZZ cited by ZZ
-- warning: missing node: SheafMap \\ SheafMap cited by \\
-- warning: missing node: SheafMap ^ ZZ cited by ^
-- warning: missing node: SheafMap ^** ZZ cited by ZZ
-- warning: missing node: SheafMap cited by Type
-- warning: missing node: SheafMap | SheafMap cited by |
-- warning: missing node: SheafMap || SheafMap cited by ||
-- warning: missing node: SheafOfRings ** SheafMap cited by **
-- warning: missing node: ZZ * SheafMap cited by *
-- warning: missing node: ZZ == SheafMap cited by Boolean
-- warning: missing node: ambient(SheafMap) cited by ambient
-- warning: missing node: coimage(SheafMap) cited by coimage
-- warning: missing node: cokernel(SheafMap) cited by cokernel
-- warning: missing node: components(SheafMap) cited by components
-- warning: missing node: cover(SheafMap) cited by cover
-- warning: missing node: degree(SheafMap) cited by ZZ
-- warning: missing node: directSum(SheafMap) cited by directSum
-- warning: missing node: dual(SheafMap) cited by dual
-- warning: missing node: dual(SheafMap,DegreeLimit=>...) cited by resolution(...,DegreeLimit=>...)
-- warning: missing node: dual(SheafMap,MinimalGenerators=>...) cited by MinimalGenerators
-- warning: missing node: dual(SheafMap,Strategy=>...) cited by Strategy
-- warning: missing node: exteriorPower(ZZ,SheafMap) cited by ZZ
-- warning: missing node: homology(SheafMap,SheafMap) cited by HH
-- warning: missing node: homomorphism'(SheafMap) cited by homomorphism'
-- warning: missing node: image(SheafMap) cited by image
-- warning: missing node: inducedMap(CoherentSheaf,CoherentSheaf,SheafMap) cited by inducedMap
-- warning: missing node: inverse(SheafMap) cited by inverse
-- warning: missing node: isIsomorphic(SheafMap,SheafMap) cited by Sequence
-- warning: missing node: isIsomorphism(SheafMap) cited by Boolean
-- warning: missing node: isIsomorphism(SheafMap) cited by isIsomorphism
-- warning: missing node: isSurjective(SheafMap) cited by Boolean
-- warning: missing node: isWellDefined(SheafMap) cited by isWellDefined
-- warning: missing node: kernel(SheafMap) cited by kernel
-- warning: missing node: lift(SheafMap) cited by lift
-- warning: missing node: map(CoherentSheaf,CoherentSheaf,SheafMap) cited by map
-- warning: missing node: matrix(SheafMap) cited by Matrix
-- warning: missing node: minimalPresentation(SheafMap) cited by minimalPresentation
-- warning: missing node: promote(SheafMap,type of Nothing) cited by promote
-- warning: missing node: promote(SheafMap,type of RingElement) cited by promote
-- warning: missing node: prune(SheafMap) cited by prune
-- warning: missing node: pullback(SheafMap,SheafMap) cited by pullback
-- warning: missing node: pushout(SheafMap,SheafMap) cited by pushout
-- warning: missing node: quotient'(SheafMap,SheafMap) cited by quotient'
-- warning: missing node: quotient(SheafMap,SheafMap) cited by quotient
-- warning: missing node: ring(SheafMap) cited by ring
-- warning: missing node: source(SheafMap) cited by source
-- warning: missing node: super(SheafMap) cited by super
-- warning: missing node: symmetricPower(ZZ,SheafMap) cited by symmetricPower
-- warning: missing node: target(SheafMap) cited by target
-- warning: missing node: tensor(SheafMap,SheafMap) cited by tensor

undocumented {
    sheafMap, -- deprecated
    -- TODO: document some of these
    (symbol SPACE, SheafMap, ZZ),
    (symbol +,   SheafMap, SheafMap),
    (symbol *,   SheafMap, SheafMap),
    (symbol *,   RingElement, SheafMap), (symbol *, ZZ, SheafMap),
    (symbol -,   SheafMap, SheafMap), (symbol -, SheafMap),
    (symbol **,  SheafMap, SheafMap),
    (symbol **,  SheafMap, SheafOfRings),  (symbol **, SheafOfRings, SheafMap),
    (symbol **,  SheafMap, CoherentSheaf), (symbol **, CoherentSheaf, SheafMap),
    (symbol |,   SheafMap, SheafMap),
    (symbol ||,  SheafMap, SheafMap),
    (symbol ==,  SheafMap, SheafMap),
    (symbol ==,  SheafMap, ZZ), (symbol ==, ZZ, SheafMap),
    }

-----------------------------------------------------------------------------
-- Types and basic constructors and methods that return a coherent sheaf map
-----------------------------------------------------------------------------

doc ///
Node
  Key
    SheafMap
  Headline
    the class of morphisms of coherent sheaves
  Description
    Text
      The most essential data of the type @TO CoherentSheaf@ in Macaulay2 is the representative
      module. For morphisms of sheaves, the data type requires a little
      more care, because even if $\mathcal{F}$ and $\mathcal{G}$ are sheaves represented by
      modules $M$ and $N$, respectively, a morphism of sheaves $\phi : \mathcal F \to \mathcal G$ is not necessarily the sheaf associated
      to a module map $\psi : M \to N$. Indeed, the best one can say is that $\phi$ is represented by some map
      $$\psi : M_{\geq d} \to N,$$
      where $d$ is some truncation degree. This means that in Macaulay2, a morphism of sheaves is
      represented as a morphism from some truncation of the source representative to the target representative.
    Example
      Q = QQ[x..z]; --Do example of accessing the truncation degree
    Text
      As illustrated in the above example, the source and target are still represented by the sheaves
      $\mathcal F$ and $\mathcal G$. The key {\tt degree} accesses the truncation degree needed to represent
      the map as a morphism of modules. To access the actual matrix representing the map, use {\tt matrix}. 
  Subnodes
    (map, CoherentSheaf, CoherentSheaf, Matrix)
    (random, CoherentSheaf, CoherentSheaf)
    (inducedMap, CoherentSheaf, CoherentSheaf)

Node
  Key
    (map, CoherentSheaf, CoherentSheaf, Matrix)
    (map, CoherentSheaf, CoherentSheaf, Matrix, InfiniteNumber)
    (map, CoherentSheaf, CoherentSheaf, Matrix, ZZ)
    (map, CoherentSheaf, Nothing,       Matrix)
    (map, Nothing,       CoherentSheaf, Matrix)
  Headline
    the constructor of morphisms of coherent sheaves

Node
  Key
    (sheaf, Matrix)
    (sheaf, Matrix, ZZ)
    (sheaf, Variety, Matrix)
    (sheaf, Variety, Matrix, ZZ)
  Headline
    the sheafification functor for morphisms

Node
  Key
    (random, CoherentSheaf, CoherentSheaf)
  Headline
    generate a random map of coherent sheaves
  Usage
    random(F, G)
  Inputs
    F:CoherentSheaf
    G:CoherentSheaf
  Outputs
    :SheafMap
--  Description
--    Text
--    Example
  SeeAlso
    setRandomSeed
    (random, Module, Module)
    (map, CoherentSheaf, CoherentSheaf, Matrix)

Node
  Key
    (inducedMap, CoherentSheaf, CoherentSheaf)
    (inducedMap, CoherentSheaf, CoherentSheaf, SheafMap)
  Headline
    induced maps on coherent sheaves

Node
  Key
    cotangentSurjection
--   (cotangentSurjection, ProjectiveVariety)

Node
  Key
    embeddedToAbstract
--   (embeddedToAbstract, ProjectiveVariety)

-----------------------------------------------------------------------------
-- Basic methods for sheaves
-----------------------------------------------------------------------------

Node
  Key
    (matrix, SheafMap)
  Headline
    the morphism of modules representing a morphisms of coherent sheaves
  Usage
    matrix phi
    phi.map
  Inputs
    phi:SheafMap
  Outputs
    :Matrix
  SeeAlso
    (truncate, List, Module)
    (truncate, Nothing, List, Matrix)

Node
  Key
    (isWellDefined, SheafMap)
  Headline
    whether a morphism of coherent sheaves is well-defined
  Usage
    isWellDefined phi
  Inputs
    phi:SheafMap
  Outputs
    :Boolean
--  Description
--    Text
--    Example
  SeeAlso
    (isWellDefined, CoherentSheaf)
    (isWellDefined, Variety)
    (isWellDefined, Matrix)

-- Note: this uses F.cache.TorsionFree and F.cache.GlobalSectionLimit,
-- which are set by twistedGlobalSectionsModule and called by HH^0(SumOfTwists) and prune(CoherentSheaf)
///

-----------------------------------------------------------------------------
-- Arithmetic operations
-----------------------------------------------------------------------------

doc ///
Node
  Key
    (quotient,  SheafMap, SheafMap)
    (quotient', SheafMap, SheafMap)
    (symbol //, SheafMap, SheafMap)
    (symbol \\, SheafMap, SheafMap)
  Headline
    factoring a morphism of coherent sheaves through another

Node
  Key
    (symbol ^, SheafMap, ZZ)
  Headline
    raises a SheafMap to a power
  Usage
    f ^ n
  Inputs
    f:SheafMap
      a sheaf map
    n:ZZ
      a non-negative integer
  Outputs
    :SheafMap
      the sheaf map f raised to the power n
  Description
    Text
      This will create a sheaf map g which is the square of f.
    Example
      R = QQ[x,y,z]
      X = variety R
      f = sheafHom(OO_X, OO_X)
      g = f ^ 2

Node
  Key
    (symbol ^**, SheafMap, ZZ)
  Headline
    repeatedly tensors a map
  Usage
    f ^** n
  Inputs
    f:SheafMap
      a sheaf map
    n:ZZ
      a non-negative integer
  Outputs
    :SheafMap
      the sheaf map f tensor n times
  Description
    Text
      This will create a sheaf map g which is the tensor of f with itself.
    Example
      R = QQ[x,y,z]
      X = variety R
      f = sheafHom(OO_X, OO_X)
      g = f ^** 2

Node
  Key
    (symbol ++, SheafMap, SheafMap)
  Headline
    takes the direct sum of two sheaf maps
  Usage
    f ++ g
  Inputs
    f:SheafMap
      a sheaf map
    g:SheafMap
      another sheaf map
  Outputs
    :SheafMap
      the direct sum of sheaf maps f and g
  Description 
    Text
        This will create a sheaf map h which is the direct sum of f and g.
    Example
      R = QQ[x,y,z]
///

///
Node
  Key
  Headline
  Usage
  Inputs
  Outputs
  Description
    Text
    Example
  SeeAlso
///
