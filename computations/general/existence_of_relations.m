load "src/StarQuotientMeasures.m";
load "src/genera_and_bounds.m";

// Main loop: verify  g_0^*(N) > B + 1,  B := sum_{p|N} B(p,N,W_N),  for all squarefree N <= N_max.
// The existence guarantee is g > B+1 (not g > B): the trace-zero normalization frees the c_0
// (identity) direction, so the thick-loop conditions live on the (g-1)-dimensional space of
// non-scalar polynomials, and at most B conditions there leave a non-scalar solution as soon as
// (g-1) - B > 0, i.e. g > B+1.  The borderline g = B+1 is inconclusive by the dimension count,
// so the exhaustive Brandt-module check below must also cover it: we run it whenever g <= B+1.
for N in [1..N_max] do
    if IsDivisibleBy(N, 10000) then
        print N;
    end if;
    if IsSquarefree(N) then
        WN := HallDivisors(N);
        g := GenusX0Nstar(N);
        if g in {0,1,2} then
            printf "N = %o: g = %o is known\n", N, g;
            continue;
        end if;
        B := BoundOnB(N, WN);
        if g le B+1 then
            // Borderline (rare): cross-check the closed-form genus against the geometric
            // Brandt-module genus before running the explicit correspondence search.
            assert StarQuotientGenus(N) eq g;
            printf "N = %o has g = %o <= %o = B+1", N, g, B+1;
            if StarQuotientHasGoodCorrespondence(N) then
                printf ", but there exists a non-trivial poly.\n";
            else
                printf ": there is no good correspondence.\n";
            end if;
        end if;
    end if;
end for;
printf "Done!\n";
