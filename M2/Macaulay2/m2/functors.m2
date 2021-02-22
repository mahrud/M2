--		Copyright 1995 by Daniel R. Grayson
-*
  functors.m2: id, HH
  ext.m2:      Ext
  tor.m2:      Tor
  varieties.m2: OO, hh
  ringmap.m2:  RingMap_*, RingMap^*
  Add:
   Tot, PP (including for line bundles?)
*-

needs "code.m2" -- for methods
needs "expressions.m2"
needs "methods.m2"

protect argument
protect subscript
protect superscript

-----------------------------------------------------------------------------
-- Local utilities
-----------------------------------------------------------------------------

wrongDomain := (G, op, X) -> error("no method for ", toString G, toString op, toString X)

-----------------------------------------------------------------------------
-- helpers for functors
-----------------------------------------------------------------------------

-- flatten the arguments given to a functor
functorArgs = method()
functorArgs(Thing,        Sequence) := (i,    args) -> prepend(i, args)
functorArgs(Thing, Thing, Sequence) := (i, j, args) -> prepend(i, prepend(j, args))
functorArgs(Thing, Thing, Thing)    :=
functorArgs(Thing, Thing)           := identity

-- check if a function can be applied to the inputs
applyFunctor = (key, G, op, X) -> (
    if (F := lookup key) =!= null then F X else wrongDomain(G, op, X))

-- TODO: retire this
applyFunctor' := (key, desc, X) -> (
    if (F := lookup key) =!= null then F X
    else error("no method for ", desc, " applied to ", X))

-----------------------------------------------------------------------------
-- Functor and ScriptedFunctor type declarations
-----------------------------------------------------------------------------

Functor = new Type of MutableHashTable
Functor.synonym = "functor"
globalAssignment Functor

ScriptedFunctor = new Type of Functor
ScriptedFunctor.synonym = "scripted functor"

-----------------------------------------------------------------------------
-- Main methods
-----------------------------------------------------------------------------
-- TODO: domain and codomain

-- TODO: should toString symbol SPACE print " " or "SPACE"?
Functor           Thing := (G, X) -> if G#?argument    then G#argument X    else wrongDomain(G, " ", X)
ScriptedFunctor _ Thing := (G, i) -> if G#?subscript   then G#subscript i   else wrongDomain(G, symbol_, i)
ScriptedFunctor ^ Thing := (G, i) -> if G#?superscript then G#superscript i else wrongDomain(G, symbol^, i)

-----------------------------------------------------------------------------
-- printing

-- TODO: improve expression and mathML of functors
net        Functor := F -> if F.?net      then F.net      else (lookup(     net, Type)) F
toString   Functor := F -> if F.?toString then F.toString else (lookup(toString, Type)) F
expression Functor := F -> new Holder from { F }
precedence Functor := x -> 70

-- TODO: use codeHelpers to get code HH to work
methods' := lookup(methods, Symbol)
methods Functor := F -> (
    if F === HH then join(methods homology, methods cohomology) else methods' F)
-- TODO: perhaps give info about argument, subscript, superscript?
methodOptions Functor := F -> null

-----------------------------------------------------------------------------
-- id: the identity morphism
-----------------------------------------------------------------------------

id = new ScriptedFunctor from {
    subscript => X -> applyFunctor((id, class X), id, symbol_, X),
    }

-----------------------------------------------------------------------------
-- OO: the structure sheaf
-----------------------------------------------------------------------------

OO = new ScriptedFunctor from {
    -- TODO: this is a weird setup :/
    subscript => X -> applyFunctor'((symbol _,     OO, class X), "OO_" | net X, (OO, X)),
    argument  => X -> applyFunctor'((symbol SPACE, OO, class X), "OO " | net X, (OO, X)),
    }
OO.texMath = ///{\mathcal O}///

-----------------------------------------------------------------------------
-- HH: homology and cohomology
-----------------------------------------------------------------------------

-- TODO: change to Options => true?
  homology = method(Options => {})
cohomology = method(Options => {Degree => 0}) -- for local cohomology and sheaf cohomology

-- TODO: is this actually necessary? functorArgs takes care of HH^i(X, Degree => ...)
  homology(ZZ, Sequence) := opts -> (i, X) ->   homology prepend(i, X)
cohomology(ZZ, Sequence) := opts -> (i, X) -> cohomology prepend(i, X)

HH = new ScriptedFunctor from {
    subscript => (
	i -> new ScriptedFunctor from {
	    -- HH_i^j X -> cohomology(i, j, X)
	    superscript => j -> new Functor from { argument => X -> cohomology functorArgs(i, j, X) },
	    -- HH_i X -> homology(i, X)
	    argument => X -> homology functorArgs(i, X)
	    }
	),
    superscript => (
	j -> new ScriptedFunctor from {
	    -- HH^j_i X -> homology(j, i, X)
	    subscript => i -> new Functor from { argument => X -> homology functorArgs(j, i, X) },
	    -- HH^j X -> cohomology(j, X)
	    argument => X -> cohomology functorArgs(j, X)
	    }
	),
    -- HH X -> homology X
    argument => X -> homology X
    }

-----------------------------------------------------------------------------
-- R^i: derived functor
-----------------------------------------------------------------------------
-- TODO

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/m2 "
-- End:
