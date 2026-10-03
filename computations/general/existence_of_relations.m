load "src/StarQuotientMeasures.m";
load "src/genera_and_bounds.m";

SetMemoryLimit(64 * 1024^3);

// A good correspondence exists if g > B+1, since at most B conditions are imposed on the
// (g-1)-dimensional space of non-scalar polynomials; for g <= B+1 we check directly.
N_max := 2589510;   // N_0 from final_inequality.m


// Class numbers for d <= N_max from the LMFDB files cl3mod8.0, cl7mod8.0, cl4mod16.0, cl8mod16.0
// (uncompressed in data/; download(ed) from
// https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups with k=0):
// h4[d] = h(-4d) for squarefree d, h1[d] = h(-d) for squarefree d = 3 mod 4, and 0 otherwise.

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

// Sanity: the known small values, the tables are nonzero exactly where they should be,
// and some random entries against ClassNumber.
assert h4[1] eq 1 and h4[5] eq 2 and h1[3] eq 1 and h1[23] eq 3;
assert h4[2] eq 1 and h4[3] eq 1 and h4[7] eq 1 and h4[11] eq 3;
for d in [1..N_max] do
    assert (h4[d] ne 0) eq IsSquarefree(d);
    assert (h1[d] ne 0) eq (d mod 4 eq 3 and IsSquarefree(d));
end for;
for i in [1..1000] do
    d := Random(1, N_max);
    if IsSquarefree(d) then
        assert h4[d] eq ClassNumber(-4*d);
        if d mod 4 eq 3 then
            assert h1[d] eq ClassNumber(-d);
        end if;
    end if;
end for;

StoreSet(_cl_store, "h4", h4);
StoreSet(_cl_store, "h1", h1);
StoreSet(_cl_store, "loaded", true);   // so _EnsureClassNumbers does not load its own tables
assert CL4(N_max) eq h4[N_max];        // CL4 reads the LMFDB table
printf "Loaded class numbers for d <= %o from the LMFDB files in data/.\n", N_max;

// A good correspondence exists if g > B+1, since at most B conditions are imposed on the
// (g-1)-dimensional space of non-scalar polynomials; for g <= B+1 we check directly.
for N in [1..N_max] do
    if IsDivisibleBy(N, 10000) then
        print N;
    end if;
    if IsSquarefree(N) then
        WN := HallDivisors(N);
        g := GenusX0Nstar(N);
        if g in {0,1} then
            printf "N = %o: g = %o is known\n", N, g;
            continue;
        end if;
        B := BoundOnB(N, WN);
        if g le B+1 then
            // Borderline (rare): cross-check the closed-form genus against the geometric
            // Brandt-module genus before running the explicit correspondence search.
            assert StarQuotientGenus(N) eq g;
            printf "N = %o: has g = %o <= %o = B+1", N, g, B+1;
            hasGoodCorrespondence, message := StarQuotientHasGoodCorrespondence(N);
            if hasGoodCorrespondence then
                printf ": there is a good correspondence: %o\n", message;
            else
                printf ": there is no good correspondence: %o\n", message;
            end if;
        end if;
    end if;
end for;
printf "Done!\n";
