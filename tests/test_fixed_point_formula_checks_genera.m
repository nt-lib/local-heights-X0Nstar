AttachSpec("local_heights.spec");

function CompareStarGenus(Nmin, Nmax)
    // Returns the squarefree N with Nmin <= N <= Nmax for which the genus of X_0(N)^*
    // computed by NTLg0star_squarefree or NTLGenusX0Nstar differs from Magma's
    // GenusX0NQuotient.
    mismatches := [];
    for N in [Nmin..Nmax] do
        if not IsSquarefree(N) then continue; end if;
        g_squarefree := NTLg0star_squarefree(N);
        g_table := NTLGenusX0Nstar(N);
        g_magma := GenusX0NQuotient(N, PrimeDivisors(N));
        if g_squarefree ne g_magma or g_table ne g_magma then
            printf "Mismatch at N = %o: NTLg0star_squarefree = %o, NTLGenusX0Nstar = %o, GenusX0NQuotient = %o\n",
                N, g_squarefree, g_table, g_magma;
            Append(~mismatches, N);
        end if;
    end for;
    return mismatches;
end function;

procedure test_star_genus_agrees_with_magma()
    TSTAssertEQ(CompareStarGenus(1, 2000), []);
end procedure;

test_star_genus_agrees_with_magma();
