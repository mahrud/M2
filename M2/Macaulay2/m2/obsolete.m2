--		Copyright 1997-2002 by Daniel R. Grayson

-- preserve this file, in case we want to remove a function and mark it obsolete:

needs "methods.m2"

-----------------------------------------------------------------------------
-- helper functions useable in documentation
-----------------------------------------------------------------------------

foo := method(Options => {})
foodict := first localDictionaries foo
---- we can get into an infinite loop by doing this: (it's like printing the contents of a mutable hash table
-- codeHelper#(functionBody value foodict#"f") = g -> {
--      ("-- method functions:", code methods value (first localDictionaries g)#"methodFunction")
--      -- ("-- option table opts:", value (first localDictionaries g)#"opts")
--      }
bar := lookup(foo,Sequence)

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/m2 "
-- End:
