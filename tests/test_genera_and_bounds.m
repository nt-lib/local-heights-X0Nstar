load "src/genera_and_bounds.m";

// The class number tables CL4(d) = h(-4d) and CL1(d) = h(-d), read from
// data/classnumbers_upto_<N_max>.dat.

procedure test_CL_known_values()
    TSTAssertEQ([CL4(d) : d in [1, 2, 3, 5, 7, 11]], [1, 1, 1, 2, 1, 3]);
    TSTAssertEQ([CL1(d) : d in [3, 23]], [1, 3]);
end procedure;

procedure test_CL_nonzero_exactly_where_defined()
    // CL4(d) is nonzero exactly for squarefree d, CL1(d) exactly for squarefree d = 3 mod 4.
    for d in [1..N_max] do
        TSTAssertEQ(CL4(d) ne 0, IsSquarefree(d));
        TSTAssertEQ(CL1(d) ne 0, d mod 4 eq 3 and IsSquarefree(d));
    end for;
end procedure;

procedure test_CL_against_ClassNumber()
    // 1000 random entries and the last one, against Magma's ClassNumber.
    for d in [Random(1, N_max) : i in [1..1000]] cat [N_max] do
        if IsSquarefree(d) then
            TSTAssertEQ(CL4(d), ClassNumber(-4*d));
            if d mod 4 eq 3 then
                TSTAssertEQ(CL1(d), ClassNumber(-d));
            end if;
        end if;
    end for;
end procedure;

test_CL_known_values();
test_CL_nonzero_exactly_where_defined();
test_CL_against_ClassNumber();
