// Run by initialize_class_numbers.sh: writes data/classnumbers_upto_<N_max>.dat, read by
// src/genera_and_bounds.m, from the LMFDB files cl3mod8.0, cl7mod8.0, cl4mod16.0, cl8mod16.0
// (uncompressed in data/; from https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups
// with k=0): h4[d] = h(-4d) for squarefree d, h1[d] = h(-d) for squarefree d = 3 mod 4,
// and 0 otherwise.
load "src/genera_and_bounds.m";   // for N_max

// Returns [h(-D) : D = r, r+m, ... <= Dmax], with 0 where -D is not fundamental, read from
// data/cl<r>mod<m>.0. Its lines "a b ..." list the fundamental -D in increasing order, with
// D = r + m*(sum of the a's so far) and b = h(-D).
function LMFDBClassNumbers(r, m, Dmax)
    file := Sprintf("data/cl%omod%o.0", r, m);
    n := (Dmax - r) div m;
    H := [Integers()| 0 : i in [0..n]];
    try
        fh := Open(file, "r");
    catch e
        error Sprintf("Cannot open %o. Download cl%omod%o.0.gz (k = 0) from https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups and uncompress it into data/.", file, r, m);
    end try;
    i := 0;
    while true do
        line := Gets(fh);
        error if IsEof(line), Sprintf("%o is incomplete: it ends before D = %o.", file, Dmax);
        entry := Split(line, "\t ");
        i +:= StringToInteger(entry[1]);
        if i gt n then break; end if;
        H[i+1] := StringToInteger(entry[2]);
    end while;
    delete fh;
    return H;
end function;

h4 := [Integers()| 0 : d in [1..N_max]];
h1 := [Integers()| 0 : d in [1..N_max]];

// d = 3 mod 4: -d is a fundamental discriminant and -4d is the discriminant of the order
// of conductor 2, so h(-4d) = (2 - (-d/2)) * h(-d) / u with u the unit index, i.e.
// h(-4d) = 3*h(-d) for d = 3 mod 8 (except h(-12) = 1, as u = 3 for d = 3) and
// h(-4d) = h(-d) for d = 7 mod 8. Reference: Cox, Theorem 7.24.
for r in [3, 7] do
    H := LMFDBClassNumbers(r, 8, N_max);
    for i in [1..#H] do
        d := r + 8*(i-1);
        h1[d] := H[i];
        h4[d] := (r eq 3) select 3*H[i] else H[i];
    end for;
end for;
h4[3] := 1;

// d = 1 mod 4 and d = 2 mod 4: -4d is a fundamental discriminant, and 4d = 4 mod 16
// resp. 4d = 8 mod 16.
for r in [4, 8] do
    H := LMFDBClassNumbers(r, 16, 4*N_max);
    for i in [1..#H] do
        h4[(r div 4) + 4*(i-1)] := H[i];
    end for;
end for;

clfile := Sprintf("data/classnumbers_upto_%o.dat", N_max);
f := Open(clfile, "w");
WriteObject(f, h4);
WriteObject(f, h1);
delete f;
printf "Wrote the class numbers for d <= %o to %o.\n", N_max, clfile;
