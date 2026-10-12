AttachSpec("local_heights.spec");
N_max := NTLNMax();

// A good correspondence exists if g > B+1, since at most B conditions are imposed on the
// (g-1)-dimensional space of non-scalar polynomials; for g <= B+1 we check directly.
for N in [1..N_max] do
    if IsDivisibleBy(N, 10000) then
        print N;
    end if;
    if IsSquarefree(N) then
        WN := NTLHallDivisors(N);
        g := NTLGenusX0Nstar(N);
        if g in {0,1} then
            printf "N = %o: g = %o is known\n", N, g;
            continue;
        end if;
        B := NTLBoundOnB(N, WN);
        if g le B+1 then
            // Borderline (rare): cross-check the closed-form genus against the geometric
            // Brandt-module genus before running the explicit correspondence search.
            assert NTLStarQuotientGenus(N) eq g;
            printf "N = %o: has g = %o <= %o = B+1", N, g, B+1;
            hasGoodCorrespondence, message := NTLStarQuotientHasGoodCorrespondence(N);
            if hasGoodCorrespondence then
                printf ": there is a good correspondence: %o\n", message;
            else
                printf ": there is no good correspondence: %o\n", message;
            end if;
        end if;
    end if;
end for;
printf "Done!\n";
