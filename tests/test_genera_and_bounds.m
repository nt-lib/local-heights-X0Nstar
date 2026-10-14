AttachSpec("local_heights.spec");
N_max := NTLNMax();

// The class number tables NTLCL4(d) = h(-4d) and NTLCL1(d) = h(-d), read from
// data/classnumbers_upto_<N_max>.dat.

procedure test_CL_known_values()
    TSTAssertEQ([NTLCL4(d) : d in [1, 2, 3, 5, 7, 11]], [1, 1, 1, 2, 1, 3]);
    TSTAssertEQ([NTLCL1(d) : d in [3, 23]], [1, 3]);
end procedure;

procedure test_CL_nonzero_exactly_where_defined()
    // NTLCL4(d) is nonzero exactly for squarefree d, NTLCL1(d) exactly for squarefree d = 3 mod 4.
    for d in [1..N_max] do
        TSTAssertEQ(NTLCL4(d) ne 0, IsSquarefree(d));
        TSTAssertEQ(NTLCL1(d) ne 0, d mod 4 eq 3 and IsSquarefree(d));
    end for;
end procedure;

procedure test_CL_against_ClassNumber()
    // 1000 random entries and the last one, against Magma's ClassNumber.
    // The d are not seeded, so a failing d is printed to make the failure reproducible.
    for d in [Random(1, N_max) : i in [1..1000]] cat [N_max] do
        if not IsSquarefree(d) then continue; end if;
        computed := [NTLCL4(d)] cat ((d mod 4 eq 3) select [NTLCL1(d)] else []);
        expected := [ClassNumber(-4*d)] cat ((d mod 4 eq 3) select [ClassNumber(-d)] else []);
        if computed ne expected then
            printf "Class number mismatch at d = %o: table %o, ClassNumber %o.\n", d, computed, expected;
        end if;
        TSTAssertEQ(computed, expected);
    end for;
end procedure;

test_CL_known_values();
test_CL_nonzero_exactly_where_defined();
test_CL_against_ClassNumber();
