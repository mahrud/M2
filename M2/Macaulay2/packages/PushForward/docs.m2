doc ///
    Key
        PushForward
    Headline
        methods to compute the pushforward of a module along a ring map
    Description
        Text
            Given a ring map $f \colon A \to B$, and a $B$-module $M$,
            $M$ has the structure of an $A$-module, and if this module is
            finitely generated over $A$, the routine @TO pushFwd@ in this package
            will compute such an $A$-module.  This is also functorial, in that if a
            map of $B$-modules (both of which are finitely generated over $A$), then
            @TO (pushFwd, RingMap, Matrix)@ will return the induced map
            on $A$-modules.

            In an algebraic sense, this is really a pull back, but thinking geometrically,
            the functions here implement the push forward of a module (or sheaf).

            This package was originally implemented by Claudiu Raicu, some changes were
            introduced by Karl Schwede, and later by David Eisenbud and Mike Stillman.
    Subnodes
        (pushFwd, RingMap)
        (pushFwd, RingMap, Module)
        (pushFwd, RingMap, Matrix)
///

doc ///
    Key
        (pushFwd, RingMap)
	(pushFwd, Ring)
    Headline
        push forward of a finite ring map
    Usage
        pushFwd f
        pushFwd B
    Inputs
        f:RingMap
	    or a ring B, and the map is taken to be the natural map from coefficientRing B
    Outputs
        :Sequence
    Description
        Text
            If $f: A \to B$ is a ring map, and $B$ is finitely generated as an $A$-module,
            then the function returns a sequence $(M, g, pf)$ containing
            (1) $M \cong B^1$ as $A$-modules,
            (2) a 1-row matrix $g$ of elements of B whose entries generate B as A-module,
            (3) a function $pf$ that
            assigns to each element $b \in B$, a matrix $A^1 \to M$,
            where the image of 1 is the element $b \in M$.
        Example
            kk = QQ;
            S = kk[a..d];
            I = monomialCurveIdeal(S, {1,3,4})
            B = S/I
            A = kk[a,d];
            f = map(B,A)
            (M,g,pf) = pushFwd f;
            M
            g
            use B
	    pf(a*b - c^2)
    Caveat
        This function is meant to be internally used.
    SeeAlso
        (pushFwd, RingMap, Module)
        (pushFwd, RingMap, Matrix)
///

doc ///
    Key
        (pushFwd, RingMap, Module)
        (pushFwd, Module)
    Headline
        push forward of a module via a finite ring map
    Usage
        N = pushFwd(f, M)
	N = pushFwd M
    Inputs
        f:RingMap
            from a ring $A$ to a ring $B$
	 	    or the natural map from coefficientRing B if f not specified
        M:Module
            a $B$-module, which via $f$ is a finite $A$-module
    Outputs
        N:Module
    Description
        Text
            Given a (not necessarily finite) ring map $f : A \to B$,
            the $B$-module $M$ can be considered as a module over $A$.
            If $M$ is finite, this method returns the corresponding
            $A$-module.
        Example
            kk = QQ;
            A = kk[t];
            B = kk[x,y]/(x*y);
            use B;
            i = ideal(x)
            f = map(B,A,{x})
            pushFwd(f,module i)
    SeeAlso
        (pushFwd, Matrix)
///

doc ///
    Key
        (pushFwd, RingMap, Matrix)
        (pushFwd, Matrix)
    Headline
        push forward of a module map via a finite ring map
    Usage
        gA = pushFwd(f, g)
        gA = pushFwd g
    Inputs
        f:RingMap
            from a ring $A$ to a ring $B$
    	 	 or (if not specified) the natural map from A = coefficientRing ring g
        g:Matrix
            (a matrix), $g : M_1 \to M_2$ of modules over $B$
    Outputs
        gA:Matrix
    Description
        Text
            If $M_1$ and $M_2$ are both finite generated as $A$-modules, via $f$, this returns the induced map
            on $A$-modules.
        Example
            kk = QQ
            A = kk[a,b]
            B = kk[z,t]
            f = map(B,A,{z^2,t^2})
            M = B^1/ideal(z^3,t^3)
            g = map(M,M,matrix{{z*t}})
            p = pushFwd(f,g)
            source p == pushFwd(f, source g)
            target p == pushFwd(f, target g)
            kerg = pushFwd(f,ker g)
            kerp = prune ker p

	    k = ZZ/32003
	    A = k[x,y]/(y^4-2*x^3*y^2-4*x^5*y+x^6-y^7)
	    A = k[x,y]/(y^3-x^7)
	    B = integralClosure(A, Keep =>{})
    	    describe B
	    f = map(B^1, B^1, matrix{{w_(3,0)}})
	    g = pushFwd(icMap A, f)
	    pushFwd(icMap A, f^2) == g*g

	    A = kk[x]
	    B = A[y, Join => false]/(y^3-x^7)
	    pushFwd B^1
	    pushFwd matrix{{y}}
        Text
	    Pushforward is linear and respects composition:
    	Example
	    B = A[y,z,Join => false]/(y^3 - x*z, z^3-y^7);
            pushFwd B^1
	    fy = pushFwd matrix{{y}}
	    fz = pushFwd matrix{{z}};
	    fx = pushFwd matrix{{x_B}};
    	    g =  pushFwd matrix{{y*z -x_B*z^2}}
	    g == fy*fz-fx*fz^2
	    fz^3-fy^7 == 0
    SeeAlso
        (pushFwd, Module)
///

document{
  Key => pushFwd,
  Headline => "push forward",
  "The push forward functor",
  UL {
       {TO (pushFwd,RingMap)," - for a finite ring map"},
       {TO (pushFwd,RingMap,Module), " - for a module"},
       {TO (pushFwd,RingMap,Matrix), " - for a map of modules"}
     }
  }

-- document{
--   Key => (pushFwd,RingMap),
--   Headline => "push forward of a finite ring map",
--   Usage => "pushFwd f",
--   Inputs => { "f" },
--   Outputs => {{"a presentation of the target of ",TT "f", " as a module over the source"},{"the matrix of generators of the target of ",TT "f"," as a module over the source"},{"a map that assigns to each element of the target of ", TT "f"," its representation as an element of the pushed forward module"}},
--   EXAMPLE lines ///
--   kk = QQ
--   S = kk[a..d]
--   I = monomialCurveIdeal(S, {1,3,4})
--   R = S/I
--   A = kk[a,d]
--   use R
--   F = map(R,A)
--   (M,g,pf) = pushFwd F;
--   M
--   g
--   pf(a*b - c^2)
--   ///,
--   Caveat => {TEX "In a previous version of this package, the third output was a function which assigned to each element of the target of ", TT "f", " its representation as an element of a free module
--       which surjected onto the pushed forward module."}
--   }

-- document{
--   Key => (pushFwd,RingMap,Module),
--   Headline => "push forward of a module",
--   Usage => "pushFwd(f,N)",
--   Inputs => { "f", "N" },
--   Outputs => {{"a presentation of ",TT "N", " as a module over the source of ",TT "f"}},
--   TEX "Given a (not necessarily finite) ring map $f:A \\rightarrow{} B$ and a $B$-module $N$ which is finite over $A$, the function returns a presentation of $N$ as an $A$-module.",
--   PARA{},
--   EXAMPLE lines ///
--   kk = QQ
--   A = kk[t]
--   B = kk[x,y]/(x*y)
--   use B
--   i = ideal(x)
--   f = map(B,A,{x})
--   pushFwd(f,module i)
--   ///
--   }

-- document{
--   Key => (pushFwd,RingMap,Matrix),
--   Headline => "push forward of a map of modules",
--   Usage => "pushFwd(f,d)",
--   Inputs => { "f", "d" },
--   Outputs => {{"the push forward of the map d"}},
--   EXAMPLE lines ///
--   kk = QQ
--   R = kk[a,b]
--   S = kk[z,t]
--   f = map(S,R,{z^2,t^2})
--   M = S^1/ideal(z^3,t^3)
--   g = map(M,M,matrix{{z*t}})
--   p = pushFwd(f,g)
--   kerg = pushFwd(f,ker g)
--   kerp = prune ker p
--   ///
--   }

doc ///
    Key
        (isModuleFinite, RingMap)
        (isModuleFinite, Ring)
        isModuleFinite
    Headline
        whether the target of a ring map is finitely generated over source
    Usage
        isModuleFinite f
        isModuleFinite R
    Inputs
        f:RingMap
            or $R$ @ofClass Ring@
    Outputs
        :Boolean
    Description
        Text
            A ring map $f \colon A \to B$ makes $B$ into a module over $A$.
            This method returns true if and only if this module is a finitely generated
            $A$-module.
        Example
            kk = QQ;
            A = kk[t];
            C = kk[x,y];
            B = C/(y^2-x^3);
            f = map(A, B, {t^2, t^3})
            isWellDefined f
            isModuleFinite f
        Example
            f = map(kk[x,y], A, {x+y})
            assert not isModuleFinite f
        Text
            If a ring $R$ is given, this method returns true if and only if $R$
            is a finitely generated module over its coefficient ring.
        Example
            A = kk[x]
            B = A[y]/(y^3+x*y+3)
            isModuleFinite B
    SeeAlso
        pushFwd
///

doc ///
Node
  Key
    [pushFwd, MinimalGenerators]
  Headline
    whether to prune the result
  Description
    Text
      This is an optional argument for the @TO pushFwd@ function. Its default value is {\tt true},
      which means that the presentation of a resulting module is @TO2(prune, "pruned")@ by default.
    Example
      R5 = QQ[a..e]
      R6 = QQ[a..f]
      M = coker genericMatrix(R6, a, 2, 3)
      G = map(R6, R5, {a+b+c+d+e+f, b, c, d, e})
      N0 = pushFwd(G, M, MinimalGenerators => false)
      N1 = pushFwd(G, M)
  SeeAlso
    prune
    trim
///
