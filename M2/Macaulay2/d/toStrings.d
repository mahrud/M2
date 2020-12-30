--		Copyright 1994 by Daniel R. Grayson
-- definitions for tostring and o << n

use arithmetic;
use stdio;

-----------------------------------------------------------------------------

putdigit(o:file, x:int):void := o << (x + if x < 10 then '0' else 'a' - 10);

putneg(o:file,x:int):void := (
    if x < 0 then (
	q := x/10;
	r := x%10;
	if r>0 then (r=r-10;q=q+1);
	putneg(o,q);
	putdigit(o,-r)));

digits(o:varstring,x:double,a:int,b:int):void := (
    x = x + 0.5 * pow(10.,double(1-a-b));
    if x >= 10. then (x = x/10.; a = a+1; b = if b==0 then 0 else b-1);
    while a > 0 do (
	putdigit(o,int(x));
	x = 10. * (x - double(int(x)));
	a = a-1;
	);
    o << '.';
    lim := pow(10.,double(-b+1));
    while b > 0 do (
	if x < lim then break;
	putdigit(o,int(x));
	x = 10. * (x - double(int(x)));
	lim = lim * 10.;
	b = b-1;
	));

-- TODO: is this really worth it??
someblanks := new array(string) len 20 at n do provide new string len n do provide ' ';

-----------------------------------------------------------------------------

export blanks(n:int):string := if n < length(someblanks) then someblanks.n else new string len n do provide ' ';

-- pads string by spaces
-- TODO: never used; delete?
export (o:file) << (s:string, n:int) : file := o << (
    if n < 0
    then ( if length(s) >= -n then s else     blanks(-n - length(s)) + s )
    else ( if length(s) >=  n then s else s + blanks( n - length(s))     ));
export (o:file) << (i:int,    n:int) : file := o << (tostring(i), n);

export     tostring(b:bool) : string := if b then "true" else "false";
export (o:file) << (b:bool) : file   := o << tostring(b);

export (o:file) << (x:int)  : file := (
    if x == 0 then putdigit(o, 0)
    else (
	if x < 0
	then ( o << '-'; putneg(o,  x); )
	else             putneg(o, -x); );
    o);
export (o:file) << (x:uchar)  : file := o << int(x);
export (o:file) << (x:short)  : file := o << int(x);
export (o:file) << (x:ushort) : file := o << int(x);
export (o:file) << (x:long)   : file := o << tostring(x);
export (o:file) << (x:ulong)  : file := o << tostring(x);

-- TODO: deprecate, use stringstream formatting
tostring5(
    x:double,						-- the number to format
    s:int,					-- number of significant digits
    l:int,					   -- max number leading zeroes
    t:int,				    -- max number extra trailing digits
    e:string			     -- separator between mantissa and exponent
    ) : string := (
    o := newvarstring(25);
    if isinf(x) then return "infinity";
    if isnan(x) then return "NotANumber";
    if x==0. then return "0.";
    if x<0. then (o << '-'; x=-x);
    oldx := x;
    i := 0;
    if x >= 1. then (
	until x < 10000000000. do ( x = x/10000000000.; i = i + 10 );
	until x < 100000.      do ( x = x/100000.;      i = i + 5 );
	until x < 100.         do ( x = x/100.;         i = i + 2 );
	until x < 10.          do ( x = x/10.;          i = i + 1 );
	)
    else (
	until x >= 1./10000000000. do ( x = x*10000000000.; i = i - 10 );
	until x >= 1./100000.      do ( x = x*100000.;      i = i - 5 );
	until x >= 1./100.         do ( x = x*100.;         i = i - 2 );
	until x >= 1.              do ( x = x*10.;          i = i - 1 );
	);
    -- should rewrite this so the format it chooses is the one that takes the
    -- least space, preferring not to use the exponent when it's a tie
    if i < 0 then (
	if -i <= l
	then   digits(o, oldx, 1, s - i - 1)
	else ( digits(o, x,    1, s - 1); o << e << tostring(i); )
	)
    else if i + 1 > s then (
	if i + 1 - s <= t
	then  digits(o, x, i + 1, 0)
	else (digits(o, x, 1, s - 1); o << e << tostring(i);)
	)
    else digits(o, x, i + 1, s - i - 1);
    tostring(o));

export   tostringRR(x:double) : string := tostring5(x,6,5,5,"e");
export (o:file) << (x:double) : file   := o << tostringRR(x);
