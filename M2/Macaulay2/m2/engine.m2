--		Copyright 1993-2002 by Daniel R. Grayson

needs "expressions.m2"
needs "integers.m2"

spliceInside = x -> new class x from deepSplice toSequence x

-- basic type

setAttribute(RawObject,ReverseDictionary,symbol RawObject)
RawObject.synonym = "raw object"
raw RawObject := x -> error "'raw' received a raw object"
net RawObject := o -> net toString o

-- monomials

setAttribute(RawMonomial,ReverseDictionary,symbol RawMonomial)
RawMonomial.synonym = "raw monomial"

RawMonomial == RawMonomial := (x,y) -> x === y
RawMonomial : RawMonomial := (x,y) -> rawColon(x,y)
ZZ == RawMonomial := (i,x) -> x == i

standardForm = method()
standardForm RawMonomial := m -> new HashTable from toList rawSparseListFormMonomial m
expression RawMonomial := x -> (
     v := rawSparseListFormMonomial x;
     if #v === 0 then expression 1
     else new Product from apply(v, (i,e) -> new Power from {vars i, e})
     )
exponents = method()
exponents(ZZ,RawMonomial) := (nvars,x) -> (
     z := new MutableList from (nvars : 0);
     scan(rawSparseListFormMonomial x, (i,e) -> z#i = z#i + e);
     toList z)
degree RawMonomial := x -> error "degree of raw monomial not defined (no monoid)"
gcd(RawMonomial,RawMonomial) := (x,y) -> rawGCD(x,y)

-- monomial orderings

setAttribute(RawMonomialOrdering,ReverseDictionary,symbol RawMonomialOrdering)
RawMonomialOrdering.synonym = "raw monomial ordering"

-- TODO: also add rawJoinMonomialOrdering
RawMonomialOrdering ** RawMonomialOrdering := RawMonomialOrdering => rawProductMonomialOrdering

-- used for debugging mgb interface, moved from ofcm.m2
monomialOrderMatrix = method()
monomialOrderMatrix RawMonomialOrdering := mo -> (
    nvars := rawNumberOfVariables mo;
    mat := rawMonomialOrderingToMatrix mo;
    -- the last entry of 'mat' determines whether the tie breaker is Lex or RevLex.
    -- there may be no other elements of mat, so the next line needs to handle that case.
    ordermat := if #mat === 3 then map(ZZ^0, ZZ^nvars, 0) else matrix pack(drop(mat, -3), nvars);
    (ordermat,
	if mat#-3 ==  0 then Lex else RevLex,
	if mat#-2 == -1 then Position => Down else
	if mat#-2 ==  1 then Position => Up   else Position => mat#-2,
	"ComponentBefore" => mat#-1)
    )

-- monoids

setAttribute(RawMonoid,ReverseDictionary,symbol RawMonoid)
RawMonoid.synonym = "raw monoid"

-- rings

setAttribute(RawRing,ReverseDictionary,symbol RawRing)
RawRing.synonym = "raw ring"
ZZ.RawRing = rawZZ()

-- ring elements (polynomials)

setAttribute(RawRingElement,ReverseDictionary,symbol RawRingElement)
RawRingElement.synonym = "raw ring element"
RawRingElement == RawRingElement := (x,y) -> x === y

RawRing _ ZZ := (R,n) -> rawRingVar(R,n)
raw Number := x -> x_((class x).RawRing)
raw InexactNumber := x -> x_((ring x).RawRing)
Number _ RawRing := (n,R) -> rawFromNumber(R,n)
RawRingElement _ RawRing := (x,R) -> rawPromote(R,x)

RawRingElement == RawRingElement := (x,y) -> x === y

ring RawRingElement := rawRing
degree RawRingElement := rawMultiDegree
denominator RawRingElement := rawDenominator
numerator RawRingElement := rawNumerator
isHomogeneous RawRingElement := rawIsHomogeneous
fraction(RawRing,RawRingElement,RawRingElement) := rawFraction

RawRingElement + ZZ := (x,y) -> x + y_(rawRing x)
ZZ + RawRingElement := (y,x) -> y_(rawRing x) + x
RawRingElement - ZZ := (x,y) -> x - y_(rawRing x)
ZZ - RawRingElement := (y,x) -> y_(rawRing x) - x
RawRingElement * ZZ := (x,y) -> x * y_(rawRing x)
ZZ * RawRingElement := (y,x) -> y_(rawRing x) * x
RawRingElement == ZZ := (x,y) -> if y === 0 then rawIsZero x else x === y_(rawRing x)
ZZ == RawRingElement := (y,x) -> if y === 0 then rawIsZero x else y_(rawRing x) === x

RawRingElement // ZZ := (x,y) -> x // y_(rawRing x)
ZZ // RawRingElement := (y,x) -> y_(rawRing x) // x

compvals := hashTable { 0 => symbol == , 1 => symbol > , -1 => symbol < }
comparison := n -> compvals#n
RawRingElement ? RawRingElement := (f,g) -> comparison rawCompare(f,g)

quotientRemainder(RawRingElement,RawRingElement) :=
quotientRemainder(ZZ, ZZ) := rawDivMod

-- monomial ideals

setAttribute(RawMonomialIdeal,ReverseDictionary,symbol RawMonomialIdeal)
RawMonomialIdeal.synonym = "raw monomial ideal"

-- free modules

setAttribute(RawFreeModule,ReverseDictionary,symbol RawFreeModule)
RawFreeModule.synonym = "raw ring"

RawFreeModule ++ RawFreeModule := rawDirectSum

degrees RawFreeModule := rawMultiDegree

ZZ _ RawFreeModule := (i,F) -> (
     if i === 0 then rawZero(F,F,0)
     else error "expected integer to be 0"
     )

RawRing ^ ZZ := (R,i) -> rawFreeModule(R,i)
RawRing ^ List := (R,i) -> (
     i = spliceInside i;
     v := - flatten i;
     if #v === 0 then rawFreeModule(R, #i ) 
     else rawFreeModule(R,toSequence v ))

rank RawFreeModule := rawRank

RawFreeModule == RawFreeModule := (v,w) -> v === w
RawFreeModule ** RawFreeModule := rawTensor

-- matrices

setAttribute(ReverseDictionary,RawMatrix,symbol RawMatrix)
RawMatrix.synonym = "raw matrix"

setAttribute(RawMutableMatrix,ReverseDictionary,symbol RawMutableMatrix)
RawMutableMatrix.synonym = "raw mutable matrix"

rawExtract = method()

rawExtract(RawMatrix,ZZ,ZZ) := 
rawExtract(RawMutableMatrix,ZZ,ZZ) := (m,r,c) -> rawMatrixEntry(m,r,c)

rawExtract(RawMatrix,Sequence,Sequence) := 
rawExtract(RawMutableMatrix,Sequence,Sequence) := (m,r,c) -> rawSubmatrix(m,spliceInside r,spliceInside c)

RawMatrix _ Sequence := 
RawMutableMatrix _ Sequence := (m,rc) -> ((r,c) -> rawExtract(m,r,c)) rc

RawMutableMatrix == RawMutableMatrix := RawMatrix == RawMatrix := rawIsEqual

RawMatrix == ZZ := RawMutableMatrix == ZZ := (v,n) -> if n === 0 then rawIsZero v else error "comparison with nonzero integer"
ZZ == RawMatrix := ZZ == RawMutableMatrix := (n,v) -> if n === 0 then rawIsZero v else error "comparison with nonzero integer"

target RawMatrix := o -> rawTarget o
source RawMatrix := o -> rawSource o
transposeSequence := t -> pack(#t, mingle t)
isHomogeneous RawMatrix := rawIsHomogeneous
entries RawMutableMatrix := entries RawMatrix := m -> table(rawNumberOfRows m,rawNumberOfColumns m,(i,j)->rawMatrixEntry(m,i,j))

ZZ * RawMatrix := (n,f) -> (
     R := rawRing rawTarget f;
     n_R * f
     )
QQ * RawMatrix := (n,f) -> (
     R := rawRing rawTarget f;
     n_R * f
     )

RawMatrix ** RawMatrix := rawTensor

rawConcatColumns = (mats) -> rawConcat toSequence mats
rawConcatRows = (mats) -> rawDual rawConcat apply(toSequence mats,rawDual)
rawConcatBlocks = (mats) -> rawDual rawConcat apply(toSequence mats, row -> rawDual rawConcat toSequence (raw \ row))

new RawMatrix from RawRingElement := (RawMatrix,f) -> rawMatrix1(rawFreeModule(ring f,1),1,1:f,0)
new RawMatrix from RawMutableMatrix := rawMatrix
new RawMutableMatrix from RawMatrix := rawMutableMatrix

RawMutableMatrix _ Sequence = (M,ij,val) -> ((i,j) -> (rawSetMatrixEntry(M,i,j,val); val)) ij

degree RawMatrix := rawMultiDegree
degrees RawMatrix :=f -> {rawMultiDegree rawTarget f,rawMultiDegree rawSource f}

clean(RR,RawMatrix) := rawClean
clean(RR,RawRingElement) := rawClean
clean(RR,RawMutableMatrix) := rawClean

norm(RR,RawMatrix) := rawNorm
norm(RR,RawRingElement) := rawNorm
norm(RR,RawMutableMatrix) := rawNorm

-- computations

setAttribute(RawComputation,ReverseDictionary,symbol RawComputation)
RawComputation.synonym = "raw computation"
status RawComputation := opts -> c -> rawStatus1 c
RawComputation_ZZ := (C,i) -> rawResolutionGetMatrix(C,i)

-- Groebner bases

RawMatrix % RawComputation := (m,g) -> rawGBMatrixRemainder(g,m)
show RawComputation := C -> rawShowComputation C

-- ring maps

setAttribute(RawRingMap,ReverseDictionary,symbol RawRingMap)
RawRingMap.synonym = "raw ring map"
RawRingMap == RawRingMap := (v,w) -> v === w

-- clean up


-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/m2 "
-- End:
