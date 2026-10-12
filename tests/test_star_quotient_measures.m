load "src/StarQuotientMeasures.m";

procedure test_dual_graph_lengths_are_integers()
    // X_0(187)^* at ell = 17 has one loop of length 2 and two of length 1 (subsec:examples_dual_graphs).
    _, stabiliser_sizes, weights := StarQuotientMeasures(17, 11, []);
    // The order of the loops depends on the Brandt module, so compare sorted lengths.
    TSTAssertEQ(Sort([stabiliser_sizes[j]*weights[j] : j in [1..#weights]]), [1, 1, 2]);
    TSTAssertEQ(Type(stabiliser_sizes[1]), RngIntElt);
end procedure;

procedure test_hecke_generator_is_checked_on_all_cusp_forms()
    // N = 37*43 (rem:Tp_generator): the newforms of 37.a1 and 43.a1 both contribute a line to
    // S_2(Gamma_0(N))^{W_N} on which T_2 acts by -2, so T_2 does not generate the Hecke algebra.
    // The old check on the dual graphs alone considered T_2 to be a generator.
    T2, T3 := Explode(StarQuotientHeckeMatrices(37*43, [2, 3]));
    TSTAssertEQ(IsHeckeGenerator(T2), false);
    TSTAssertEQ(IsHeckeGenerator(T3), true);
end procedure;

procedure test_relations_only_for_visible_loops()
    // X_0(187)^* has a single visible loop: the loop of length 2 at ell = 17. At ell = 11 all
    // loops have length 1, which only contributes a zero row.
    relations, is_hecke_generator, hecke_on_S := StarQuotientHeightRelations(187, [2]);
    TSTAssertEQ(#[row : row in relations[1] | row ne [0, 0, 0]], 1);
    TSTAssertEQ(is_hecke_generator, [true]);
    TSTAssertEQ(Nrows(hecke_on_S[1]), 3);
end procedure;

procedure test_has_good_correspondence()
    // Levels from sec:examples, and N = 85, where there is no good correspondence.
    TSTAssertEQ([StarQuotientHasGoodCorrespondence(N) : N in [187, 247, 319, 85]],
                [true, true, true, false]);
end procedure;

procedure test_has_good_correspondence_does_not_hide_errors()
    // Errors other than "no Hecke generator" must propagate (here: N is not squarefree),
    // instead of being swallowed while searching for a Hecke generator.
    TSTAssertRaises(StarQuotientHasGoodCorrespondence, "Assertion failed", 12);
end procedure;

test_dual_graph_lengths_are_integers();
test_hecke_generator_is_checked_on_all_cusp_forms();
test_relations_only_for_visible_loops();
test_has_good_correspondence();
test_has_good_correspondence_does_not_hide_errors();
