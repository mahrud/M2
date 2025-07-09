-- Documentation for the package name key

doc ///
    Key
        CpMackeyFunctors
    Headline
        for Cp-Mackey functor computations
    Description
        Text
            This {\em Macaulay2} package implements basic homological algebra for {\em Mackey functors} of a cyclic group $C_p$ of prime order. Some highlights of our functionality include:

            @UL{
                ("free ", TO2((resolution,CpMackeyFunctor,ZZ),"resolutions"), " of ", TEX"$C_p$","-Mackey functors"),
                (TO2((Ext,ZZ,CpMackeyFunctor,CpMackeyFunctor),"Ext"), " and ", TO2((Tor,ZZ,CpMackeyFunctor,CpMackeyFunctor),"Tor"), " computations for (cohomological) Mackey functors"),
                ("the ", TO2(boxProduct,"box product"), " and ", TO2(internalHom,"internal hom")),
                (TO2("constructing examples of Mackey functors", "constructions")," of various common examples of Mackey functors"),
                ("the ability to generate ", TO2(makeRandomCpMackeyFunctor,"random Mackey functors")," over ", TEX"$C_p$"),
                (TO2("explicit applications of the CpMackeyFunctors package","explicit examples")," from the equivariant algebra literature")
            }@

            {\bf NOTE} This package uses some methods defined in the new Complexes package. You may see warnings when you load this package, e.g.
        Pre
            warning: symbol "GradedModuleMap" in OldChainComplexes.Dictionary is shadowed by a symbol in Complexes.Dictionary
               use the synonym OldChainComplexes$GradedModuleMap
            warning: symbol "GradedModule" in OldChainComplexes.Dictionary is shadowed by a symbol in Complexes.Dictionary
               use the synonym OldChainComplexes$GradedModule
            warning: symbol "res" in OldChainComplexes.Dictionary is shadowed by a symbol in Complexes.Dictionary
               use one of the synonyms OldChainComplexes$res, OldChainComplexes$resolution
            warning: symbol "resolution" in OldChainComplexes.Dictionary is shadowed by a symbol in Complexes.Dictionary
               use one of the synonyms OldChainComplexes$res, OldChainComplexes$resolution
        Text
            This is nothing to worry about! To remove these warnings, you can add {\tt HomologicalAlgebraPackage = "Complexes"} to your initialization file. See @HREF{"https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/OldChainComplexes/html/index.html","here"}@ for more information.

    SeeAlso
        "background on Mackey functors"
        "the abelian category of Mackey functors"
        "constructing examples of Mackey functors"
        "list of common Mackey functors"
        "explicit applications of the CpMackeyFunctors package"
///