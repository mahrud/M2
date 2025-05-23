--this file is in the public domain

newPackage(
    "SectionRing",
    Date => "May 22, 2025",
    Version => "0.3",
    Headline => "the section ring of a Weil Divisor",
    Authors => {
	{ Name => "Andrew Bydlon",  Email => "thelongdivider@gmail.com", HomePage => "http://www.math.utah.edu/~bydlon/" },
	{ Name => "Mahrud Sayrafi", Email => "mahrud@umn.edu",           HomePage => "https://math.umn.edu/~mahrud" }
    },
    PackageExports => {
	"Varieties",
	"WeilDivisors",
    },
    Keywords => { "Commutative Algebra", "Projective Algebraic Geometry" },
    DebuggingMode => true,
)

export{
	"globallyGenerated",
	"isMRegular",
	"mRegularity",
	"mRegular" => "mRegularity",
	"sectionRing",
}

-----------------------------------------------------------------------
-- Find m such that OO_X(mD) is a globally generated line bundle
-----------------------------------------------------------------------

globallyGenerated = method()
globallyGenerated Ideal := I -> globallyGenerated divisor I
globallyGenerated WeilDivisor := D -> (
    -- compute the smallest positive number (using a binary search)
    -- such that OO_X(mD) is globally generated for an ample OO_X(D).
	a:=1;

	while ((1%(baseLocus(a*D)) == 0) != true) do (
		a =2*a;
	);

	upperbound := a;
	lowerbound := ceiling(a/2);

	while (lowerbound < upperbound-1) do (
		a = ceiling((lowerbound + upperbound)/2);
		if ((1%(baseLocus(a*D)) == 0) != true) then (
			lowerbound = a;
		)
		else if ((1%(baseLocus(a*D)) == 0) == true) then (
			upperbound = a;
		);

	);
	upperbound
)

-----------------------------------------------------------------------
-- Castelnuovo-Mumford's m-regularity with respect to an ample bundle B
-----------------------------------------------------------------------

isMRegular = method()
isMRegular(CoherentSheaf, ZZ) := (F, m) -> (
    -- whether a sheaf F on X is m-regular relative to OO_X(1)
    -- TODO: use Tate resolutions
    all(1 .. dim variety F, i -> HH^i(F(m-i)) == 0)
)
isMRegular(CoherentSheaf, CoherentSheaf, ZZ) := (F, B, m) -> (
    -- whether a coherent sheaf F on a projective variety X is m-regular
    -- with respect to a globally generated ample line bundle B
    -- see Definition 1.8.4 in Positivity in Algebraic Geometry I
    n := dim variety F;
    G := F ** B ^** (m - n - 1);
    all(n, j -> HH^(n-j)(G **= B) == 0)
)

-----------------------------------------------------------------------

mRegularity = method()
mRegularity Ideal := I -> (
    -- compute m for which OO_X(D) is m-regular relative to OO_X(1),
    -- where D is the divisor corresponding to a codim 1 ideal I.
    mRegularity dual sheaf I
)
mRegularity CoherentSheaf := F -> (
    -- compute m for which a sheaf F in m-regular relative OO_X(1)
    V := variety F;
    mRegularity(F, OO_V(1))
)
mRegularity(CoherentSheaf, CoherentSheaf) := (F, G) -> (
    -- computes m for which a sheaf F is m-regular relative to G,
    -- in the sense of Castelnuovo-Mumford, using a binary search
	bool0 := isMRegular(F,G,0);
	m:=0;
	lowerbound:=0;
	upperbound:=0;
	a:=0;

	if (bool0 == true) then (
--Tests for a negative-regularity in the case that F is 0-regular relative to G
		m=-1;
		while (isMRegular(F,G,m)) do (
			m=2*m;
		);

		lowerbound = m;
		upperbound = ceiling(m/2);

		while (lowerbound < upperbound-1) do (
			a = ceiling((lowerbound + upperbound)/2);
			if (isMRegular(F,G,a) != true) then (
				lowerbound = a;
			)
			else if (isMRegular(F,G,a) == true) then (
				upperbound = a;
			);
		);
	)

	else if (bool0==false) then (
--Tests for positive-regularity in the case that F is NOT 0-regular relative to G
		m=1;
		while (isMRegular(F,G,m) != true) do (
			m=2*m;
		);

		upperbound = m;
		lowerbound = ceiling(m/2);

		while (lowerbound < upperbound-1) do (
			a = ceiling((lowerbound + upperbound)/2);
			if (isMRegular(F,G,a) != true) then (
				lowerbound = a;
			)
			else if (isMRegular(F,G,a) == true) then (
				upperbound = a;
			);
		);
	);
	upperbound
)

-----------------------------------------------------------------------
-- Ring of sections of an ample line bundle on a projective variety
-----------------------------------------------------------------------

-- two small utilities for working with vectors of scalars
-- TODO: there must be a better way to implement these
substituteScalarVector = (R, L) -> apply(L, z -> sub(z, R))
isScalarVector = L -> (R := ring(L#0); all(L, z -> z == 0 or degree z == degree 1_R))

-- Given an ideal $I$ as input, dualizes the ideal, and maps it back into the ring,
-- producing $\operatorname{Hom}_R(I,R) \cong J \subset R$.
-- This method is used to produce the global sections $H^0(mD)$,
-- where $D$ is an integral divisor defined by $I$.
dualToIdeal = method()
dualToIdeal Ideal := I -> embedAsIdeal(dual module I, IsGraded => true)

sectionRing = method()
sectionRing WeilDivisor := D -> sectionRing ideal D
sectionRing Ideal := I -> (
    -- compute the ring of sections of a semi-ample divisor associated to I

	local L;
	local Rel;
	local KerT;
	local Part;
	local LengP;
	local LengPa;
	local numDegs;
	local AdmPart;
	local NumCols;
	local b;
	local e;

	R := ring(I);
					
--To Apply the Regularity Theorem of Mumford, the sheaf needs to be Globally Generated Sheaf. Thus in the case O_X(D) is not globally generated, we consider F=O_X(D),O_X(2D), ... , O_X((l-1)D) (which correspond to J#1, J#2,...) and F being relatively G-m-regular, where G = O_X(lD) is globally generated.  This produces bound, where all generators are found in lower degrees than bound.

	l := globallyGenerated(I);
	bound := l;
	G := first entries gens I;
	J :={0};

	j:=1;
	while(j<l+1) do (						
		J = J|{Hom(ideal(apply(G, z->z^j)),R)};
		j=j+1;
	);
	
	j=1;
	while(j<(l+1)) do (	
		bound = max(bound,l*(mRegularity(sheaf(J#j),sheaf(J#l)))+j);
		j = j+1;
	);
	bound = bound + 1;

--The next block of code produces a polynomial ring S with generators in degrees 1,2,3,...,bound which will then be quotiented to produce the section ring.  Here Map_i Represents the map H^0(O_X(iD)) -> J^(i) and n_i is the rank of H^0(O_X(iD)).

	KK:= coefficientRing(R);
	Z := dualToIdeal(I);
	Shift := (Z#1)#0;
	J = {0,reflexify((Z#0))};
	FF := {0,basis(Shift,J#1)};
	n := {0,numColumns(FF#1)};
	F := {0,map(R^(numRows(FF#1)),R^(n#1),FF#1)};
	Map := {0,(gens J#1)*(F#1)};
	Y := local Y;
	myVars := {toList(Y_{1,1}..Y_{1,n#1})};						
	DegreeList :={};
	l=0;
	while(l<n#1) do(
		DegreeList = DegreeList | toList({1});
		l = l+1;
	);					
	i:=2;
	while (i < bound) do(
		J = J | {reflexivePower(i,J#1)};
		FF = FF | {basis((Shift*i),J#i)};
		n = n | {numColumns((FF#i))};			
		F = F | {map(R^(numRows((FF#i))),R^(n#i),FF#i)};
		Map = Map | {(gens J#i)*(F#i)};
		myVars = myVars | {toList(Y_{i,1}..Y_{i,n#i})};  
		l=0;
		while(l<n#i) do(
			DegreeList = DegreeList | toList({i});
			l = l+1;
		);
		i=i+1;
	);

	Vars := flatten myVars;
	numVars:= #Vars;

	S := KK [Vars,Degrees=>DegreeList];
	myVars = apply(myVars, z->apply(z,x->value(x)));
	numDegs = #myVars;

--The following block of code is used to compute the relations on S which define the section ring SR.  It does so by going degree by degree (starting at degree 2), and considering the morphisms \oplus_{i=0,\ldots,[j/2]} H^0((j-i)D) \otimes H^0(iD) --> H^0(jD), computing its kernel, and multiplying the matrix representing the kernel with the corresponding vector of variables of S, Vect_{j-i} \otimes Vect_i.  This gives relations, then inserted into RelIdeal.

	RelIdeal := ideal(0);
	Spar := S;
	Vect := {0};

	c:=1;
	while((c<bound) and (n#c>0)) do (
		Vect = Vect | {transpose matrix{myVars#(c-1)}};
		c=c+1;
	);

	j=2;
	while ( (dim(Spar) >  dim(R)) or (isDomain(Spar) != true)) do (	

--Relations are achieved by finding the kernels of the direct sums of tensor products of global sections.  However, some efficiency improvements can be achieve by considering a minimal number of such sums/tensors.  To do this, I consider partitions of the given degree of interest and throw out any partitions which are either above the bound in which our generators are considered, or can be factored through another partition.  For example, O_V(D)^{\otimes 3} -> O_V(2D)\otimes O_V(D) -> O_V(3D), so if bound>1, then the partition (1,1,1) of 3 is excluded.  Throughout, a is an index for which partition is chosen, and b is an index for an element of a given partition.  Additionally, MapTot is the total map of lower degree tensors into the degree in which relations are being considered, VectTot the corresponding vector of variables.

		Part = partitions(j);
		LengP = #(Part);
		a:=0;
		AdmPart = {};
		while (a<LengP) do(
			if((Part#a#0 < bound) and ((Part#a)#(#(Part#a)-1) + (Part#a)#(#(Part#a)-2) > min(bound,j)-1)) then (
				AdmPart = AdmPart | {(Part)#a}; 
			);
			a=a+1;
		);
	
		LengP = #(AdmPart);

		a=0;

		TotMapTemp := Map#(AdmPart#a#0);
		TotVectTemp := Vect#(AdmPart#a#0);
		b=1;
		LengPa = #((AdmPart)#a);
		while (b<LengPa) do (
			TotMapTemp = TotMapTemp ** Map#(AdmPart#a#b);
			TotVectTemp = TotVectTemp ** Vect#(AdmPart#a#b);
			b = b+1;
		);

		MapTot := TotMapTemp;
		VectTot := TotVectTemp;

		a=1;
		while(a<LengP) do ( 
			TotMapTemp = Map#(AdmPart#a#0);
			TotVectTemp = Vect#(AdmPart#a#0);

			b=1;
			LengPa = #(AdmPart#a);
			while (b<LengPa) do (
				TotMapTemp = TotMapTemp ** Map#(AdmPart#a#b);
				TotVectTemp = TotVectTemp ** Vect#(AdmPart#a#b);
				b = b+1;
			);

			MapTot = MapTot | TotMapTemp;
			VectTot = VectTot || TotVectTemp;
			a = a+1;
		);

		KerT = generators ker(MapTot);		

		NumCols = numColumns(KerT);
		e = 0;
			
		while (e < NumCols) do (
			L = flatten entries KerT_{e};
			if isScalarVector L then (
				L = substituteScalarVector(S,L);
				Rel = sub((entries (matrix{L}*VectTot))#0#0,S);
				RelIdeal = trim(RelIdeal + ideal(Rel));
				Spar = S/RelIdeal;
			);
			e=e+1;
		); 
		j=j+1;
	);

--Some code to improve the presentation of the ring, both in terms of having a more standard list of generators A_1...A_N, and eliminating redundant generators

	A := local A;
	BetterS := KK[A_1..A_numVars,Degrees=>DegreeList];
	BetterMap := map(BetterS,S,toList(A_1..A_numVars));	
	BetterRelIdeal := BetterMap(sub(RelIdeal,S));
	minimalPresentation(BetterS/BetterRelIdeal)
)

-----------------------------------------------------------------------

beginDocumentation();

doc ///
Node
  Key
    SectionRing
  Headline
    computing the section ring of a Weil Divisor
  Description
    Text
      This package provides algorithms for computing the ring of sections a semi-ample Weil divisor.
    Tree
      :Main algorithm
       > sectionRing
      :Positivity Computations
       > globallyGenerated
       > mRegularity
       > isMRegular
  SeeAlso
    "WeilDivisors :: WeilDivisors"

Node
  Key
    globallyGenerated
   (globallyGenerated, WeilDivisor)
   (globallyGenerated, Ideal)
  Headline
    find smallest integer a such that OO_X(mD) is globally generated
  Usage
    globallyGenerated(D)
  Inputs
    D:{WeilDivisor,Ideal}
  Outputs
    :Number
  Description
    Text
      This method uses a binary search to find the smallest integer $m$ with the property
      that $|mD|$ is a basepoint-free linear series. In this case, the corresponding line
      bundle is globally generated.

Node
  Key
    isMRegular
   (isMRegular, CoherentSheaf, ZZ)
   (isMRegular, CoherentSheaf, CoherentSheaf, ZZ)
  Headline
    whether F is m-regular in the sense of Castelnuovo-Mumford
  Usage
    isMRegular(F,m)
    isMRegular(F,B,m)
  Inputs
    F:CoherentSheaf
      over a projective variety $X$.
    B:CoherentSheaf
      which is a globally generated ample line bundle on $X$;
      if omitted, assumes $B = \mathcal{O}_X(1)$.
    m:ZZ
  Outputs
    :Boolean
      whether $\mathcal F$ is $m$-regular with respect to $B$ in the sense of Castelnuovo-Mumford
  Description
    Text
      This method tests whether
      $$ H^i(\mathcal F \otimes B^{\otimes(m-i)}) = 0 \quad \text{for every} \quad i > 0.$$
      As soon as a non-zero cohomology is found, the algorithm stops and returns false.
      If none is found, $\mathcal F$ is $m$-$B$-regular, and it returns true.
  References
    See Definition 1.8.4 of Lazarsfeld's Positivity in Algebraic Geometry I.
  Caveat
    In the case $B = \mathcal{O}_X(1)$, it may be faster to use @TO "BGG::BGG"@
    or @TO "TateOnProducts::TateOnProducts"@ to compute many cohomologies at once.

Node
  Key
    mRegularity
   (mRegularity, Ideal)
   (mRegularity, CoherentSheaf)
   (mRegularity, CoherentSheaf, CoherentSheaf)
  Headline
    compute the Castelnuovo-Mumford regularity of F with respect to G
  Usage
    mRegularity(F)
    mRegularity(F,G)
  Inputs
    F:{CoherentSheaf,Ideal}
      over a projective variety $X$;
      given an ideal, assumes $OO_X(D)$ where $D$ is divisor associated to $I$.
    B:CoherentSheaf
      which is a globally generated ample line bundle on $X$;
      if omitted, assumes $B = \mathcal{O}_X(1)$.
  Outputs
    :ZZ
  Description
    Text
      This method utilizes a binary search to compute the smallest $m$ such that $\mathcal F$
      is $m$-regular with respect to $B$, utilizing the function @TO "isMRegular"@.
      computes the regularity of O_X(D), where D is the associated divisor to I.

Node
  Key
    sectionRing
   (sectionRing, WeilDivisor)
   (sectionRing, Ideal)
  Headline
    compute the section ring of an ample divisor
  Usage
    sectionRing(D)
  Inputs
    D:{WeilDivisor,Ideal}
      an ample divisor on a projective variety $X$;
      given an ideal, the corresponding divisor is used.
  Outputs
    :Ring
  Description
    Text
      This algorithm begins by computing the regularity $m$ of
      $\mathcal O_X, \mathcal O_X(D), \mathcal O_X(2D), \dots, \mathcal O_X((l-1)D)$
      with respect to $\mathcal O_X(lD)$, where $l$ is the output of @TT "globallyGenerated(D)"@.
      Mumford's Theorem (1.8.5 in Positivity in Algebraic Geometry I) implies that each of the maps
      $\mathcal O_X(iD)\otimes \mathcal O_X(lD)^\otimes m \to \mathcal O_X((i+ml)D)$ is surjective.
      Thus, all generators for the section ring can be assumed in lower degree than bound.
      Thus it forms a polynomial ring $S$ over the base field with $h^0(iD)$-many generators in degree $i$,
      for $i = 1,2,\dots,\mathrm{bound}-1$.

      Next, relations in degree $d$ are computed by considering the total maps
      $\oplus_{\mathrm{partitions } P \mathrm{ of } d} \otimes_{i\in P} \mathcal O_X(i D) \to \mathcal O_X(dD)$.
      Each of these relations is then quotiented, until the point that a domain of the correct dimension is produced.
      Some steps are then performed to make the output more readable and standard.
  SeeAlso
    globallyGenerated
///

-----------------------------------------------------------------------

TEST ///
R = QQ[x,y,z]/ideal(x^3+y^3-z^3);
I = ideal(x,y-z);
assert( globallyGenerated(I) == 2)
///


TEST ///
R = QQ[x,y,z]/ideal(x^4+y^4-z^4);
I = ideal(x,y-z);
assert( globallyGenerated(I) == 3)
///

TEST ///
R = QQ[x,y,z]/ideal(x^5+y^5-z^5);
I = ideal(x,y-z);
assert( globallyGenerated(I) == 4)
///

TEST ///
X = Proj(QQ[x,y,z,w,f]);
F = OO_X(4);
G = OO_X(-2);
assert(not isMRegular(F, -5) and isMRegular(F, -4))
assert(not isMRegular(F, OO_X(1), -5) and isMRegular(F, OO_X(1), -4))
assert(mRegularity F == -4 and mRegularity G == 2)
///

TEST ///
R = QQ[x,y,z,w,f];
I = ideal(x-y+w);
S = sectionRing I;
assert( (#(vars S)) == dim S)
///

TEST ///
R = QQ[x,y,z]/(x^3+y^3-z^3)
I = ideal(x,y-z);
S = sectionRing I;
J = ideal(S);
L = first entries gens J;
assert((#L==1) and ((degree(L#0))#0 == 6))
///

end

restart
needsPackage "SectionRing"
check "SectionRing"

installPackage "SectionRing"
viewHelp oo
