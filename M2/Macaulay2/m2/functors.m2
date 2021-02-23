--		Copyright 1995 by Daniel R. Grayson

needs "code.m2" -- for methods
needs "expressions.m2"
needs "methods.m2"

-----------------------------------------------------------------------------
-- helpers for functors
-----------------------------------------------------------------------------

-- flatten the arguments given to a scripted functor
functorArgs = method()
functorArgs(Thing,        Sequence) := (i,    args) -> prepend(i, args)
functorArgs(Thing, Thing, Sequence) := (i, j, args) -> prepend(i, prepend(j, args))
functorArgs(Thing, Thing, Thing)    :=
functorArgs(Thing, Thing)           := identity

wrongDomain := (G, op, X) -> error("no method for ", toString G, toString op, toString X)

-- check if a function can be applied to the inputs
applyMethod = (key, G, op, X) -> (
    if (F := lookup key) =!= null then F X else wrongDomain(G, op, X))

-- TODO: combine these with applyMethod and retire these
applyMethod' = (key, desc, X) -> (
    if (F := lookup key) =!= null then F X
    else error("no method for ", desc, " applied to ", X))

-- used in hh
applyMethod'' = (F, X) -> (
    -- TODO: write a variation of lookup to do this
    key := prepend(F, delete(Option, apply(X, class)));
    applyMethod'(key, toString F, if #X == 1 then X#0 else X))

-- TODO: combine for all functors
-- used in Ext
applyMethodWithOpts' = (key, desc, X, opts) -> (
    if (F := lookup key) =!= null then (F opts) X
    else error("no method for ", desc, " applied to ", X))

applyMethodWithOpts'' = (F, X, opts) -> (
    -- TODO: write a variation of lookup to do this
    key := prepend(F, apply(X, class));
    applyMethodWithOpts'(key, toString F, X, opts))

-----------------------------------------------------------------------------
-- Functor and ScriptedFunctor type declarations
-----------------------------------------------------------------------------

protect argument
protect subscript
protect superscript

Functor = new Type of MutableHashTable
Functor.synonym = "functor"
globalAssignment Functor

ScriptedFunctor = new Type of Functor
ScriptedFunctor.synonym = "scripted functor"

-----------------------------------------------------------------------------
-- Main methods
-----------------------------------------------------------------------------
-- TODO: domain and codomain

Functor           Thing := (G, X) -> if G#?argument    then G#argument X    else wrongDomain(G, " ", X)
ScriptedFunctor _ Thing := (G, i) -> if G#?subscript   then G#subscript i   else wrongDomain(G, symbol _, i)
ScriptedFunctor ^ Thing := (G, i) -> if G#?superscript then G#superscript i else wrongDomain(G, symbol ^, i)

-----------------------------------------------------------------------------
-- printing and introspection

expression Functor := F -> if F.?expression then F.expression else (lookup(expression, Type)) F
describe   Functor := F -> if F.?describe   then F.describe   else (lookup(describe,   Type)) F
net        Functor := F -> if F.?net        then F.net        else (lookup(net,        Type)) F
toString   Functor := F -> if F.?toString   then F.toString   else (lookup(toString,   Type)) F
toExternalString Functor := toString @@ describe
precedence Functor := x -> 70

-- TODO: use codeHelpers to get code HH to work
-- TODO: get methods OO to work
-- TODO: improve this for RingMap_*
methods' := lookup(methods, Symbol)
methods Functor := F -> (
    if F === HH then join(methods homology, methods cohomology) else methods' F)
-- TODO: perhaps give info about argument, subscript, superscript?
methodOptions Functor := F -> null

-----------------------------------------------------------------------------
-- id: the identity morphism
-----------------------------------------------------------------------------

id = new ScriptedFunctor from {
    subscript => X -> applyMethod((id, class X), id, symbol_, X),
    }

-----------------------------------------------------------------------------
-- HH: homology and cohomology
-----------------------------------------------------------------------------

-- TODO: change to Options => true?
  homology = method(Options => {})
cohomology = method(Options => {Degree => 0}) -- for local cohomology and sheaf cohomology

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

-- TODO: is this actually necessary? functorArgs takes care of HH^i(X, Degree => ...)
  homology(ZZ,Sequence) := opts -> (i,X) ->   homology prepend(i,X)
cohomology(ZZ,Sequence) := opts -> (i,X) -> cohomology(prepend(i,X), opts)

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/m2 "
-- End:
