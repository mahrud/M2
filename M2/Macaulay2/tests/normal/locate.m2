f = () -> if ( true ) then ( "true" ) else ( "false" )
assert( ((F = functionBody f;)) === null );
assert( ((c = pseudocode   f;)) === null );
assert( ((C = pseudocode   F;)) === null );
getcols = P -> (toList P)_{2,4,6}
assert( (getcols locate f) === {4,54,7} );
assert( (getcols locate f) === {4,54,7} );
assert( (getcols locate F) === {4,54,7} );
assert( (getcols locate c) === {10,54,43} );
assert( (getcols locate C) === {10,54,43} );
assert( (getcols locate C_0) === {15,19,15} );
assert( (getcols locate C_1) === {29,35,29} );
assert( (getcols locate C_2) === {45,52,45} );
f = () -> A B C
assert( ((C = pseudocode functionBody f;)) === null );
assert( (getcols locate f) === {4,15,7} );
assert( (getcols locate C) === {10,15,11} );
assert( (getcols locate C_0) === {10,11,10} );
assert( (getcols locate C_1) === {12,15,13} );
assert( (getcols locate C_1_0) === {12,13,12} );
assert( (getcols locate C_1_1) === {14,15,14} );
f = () -> (A B) C
assert( ((C = pseudocode functionBody f;)) === null );
assert( (getcols locate f) === {4,17,7} );
assert( (getcols locate C) === {10,17,11} );
assert( (getcols locate C_0) === {11,14,12} );
assert( (getcols locate C_1) === {16,17,16} );
assert( (getcols locate C_0_0) === {11,12,11} );
assert( (getcols locate C_0_1) === {13,14,13} );
end--
-*
-- to update this file simply run these lines:
src = last separate("end-{2,}", get "locate.m2");
"locate.m2" << generateAssertions src << endl;
"locate.m2" << "end" << "--" << src << close;
*-
--  0  0  1  1 1  1 2 2    2 2    3 3 3    4 4     5 5
--  4  7  0  3 5  8 0 2    7 9    4 6 8    3 5     1 3
f = () -> if ( true ) then ( "true" ) else ( "false" )
(F = functionBody f;)
(c = pseudocode   f;)
(C = pseudocode   F;)
getcols = P -> (toList P)_{2,4,6}
getcols locate f
getcols locate f-* why is focus on '->'? *-
getcols locate F
getcols locate c-* why is focus on '( "false" )'? *-
getcols locate C
getcols locate C_0
getcols locate C_1
getcols locate C_2

-- Adjacent
--        1 1 1
--  4  7  0 2 4
f = () -> A B C
(C = pseudocode functionBody f;)
getcols locate f
getcols locate C
getcols locate C_0
getcols locate C_1
getcols locate C_1_0
getcols locate C_1_1
--        1  1  1
--  4  7  0  3  5
f = () -> (A B) C
(C = pseudocode functionBody f;)
getcols locate f
getcols locate C
getcols locate C_0
getcols locate C_1
getcols locate C_0_0
getcols locate C_0_1

-- TODO: add tests for other kinds of code as well
