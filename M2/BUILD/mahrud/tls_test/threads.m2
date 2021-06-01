debug Core

N = 0
threadLocal myN'
myN' = 0


for i to 25 do (
    schedule((c) -> (
	    if myN' === null then myN' = 0;
	      N =   N + 1;
	    myN'= myN'+ 1;
	    s := "Value for " | c | ":\t" | toString N | "\tvs " | toString myN';
	    myprint s -- toString (s | "\t" | net random(ZZ^2,ZZ^2))
	    ),
	ascii(i + first ascii "A")
	);
    );

nanosleep(10000000 * random(10));

  N =   N + 1;
myN'= myN'+ 1;
s := "Value for main: " | toString N | "\tvs " | toString myN';
myprint s -- toString (s | " " | net random(ZZ^2,ZZ^2))

exit
