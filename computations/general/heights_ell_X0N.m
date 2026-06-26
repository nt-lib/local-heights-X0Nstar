N := 3 * 61;
assert IsSquarefree(N);

p := 5; // we take the Hecke correspondence T_p

printf "N = %o = %o\n", N, Factorization(N);
printf "values of local heigths h_ell for X_0(N):\n";

S := CuspidalSubspace(ModularSymbols(N, 2, +1));
trTp := Trace(HeckeOperator(S, p)) / Dimension(S);

for ell in PrimeDivisors(N) do 
    printf "ell = %o:\n", ell;
    Nprime := ExactQuotient(N, ell);

    // we take the correspondence F = T_p - tr(T_p) * I

    // compute the values on the two components (vertices) X_0(N/ell) otimes F_ell
    Sold := CuspidalSubspace(ModularSymbols(Nprime, 2, +1));
    v1 := Trace(HeckeOperator(Sold, p)) - Dimension(Sold) * trTp; // iota_1
    v2 := /*Trace(HeckeOperator(Sold, p))*/ - Dimension(Sold) * trTp; // iota_ell
    vertices := [v1, v2];
    printf "vertices: %o\n", vertices;

    // compute the values on the edges
    Smod := SupersingularModule(ell, Nprime);
    Tp := HeckeOperator(Smod, p);
    Tp_minus_ap := Tp - Identity(Parent(Tp)) * trTp;
    F := ChangeRing(Tp_minus_ap, Rationals());
    print F;

    // compute mu_F for a multiple banana graph, see [Betts--Dogra, 12.2.3]
    n := Nrows(F);
    length := func<e | MonodromyWeights(Smod)[e]>;
    nu := &+[length(e) : e in [1..n]]^-1;
    e_vector := func<e | Vector(Rationals(), [f eq e select 1 else 0 : f in [1..n]])>;
    edges := [1/length(e) * ((e_vector(e) - &+[1/(nu * length(f)) * e_vector(f) : f in [1..n]]) * F)[e] : e in [1..n]];
    printf "edges: %o\n", edges;

    wps := [AtkinLehnerOperator(Smod, d) : d in PrimeDivisors(N)];
    wds := [&*[wps[i] : i in I] : I in Subsets({1..#wps}) | #I ge 1];
    printf "%o out of %o edges are not fixed by any w_prime.\n", #{e : e in [1..n] | forall{wd : wd in wds | wd[e,e] ne 1}}, n;

    // compute j_Gamma by inverting the Laplacian on banana graph
    // j_Gamma on edge e is g_e(x) = a_e x^2 + b_e x + c_e
    // -2a_e = -g_e'' = edges[e]
    // on v_1, v_2: sum_e (g_e'(0) - g_e'(1)) = sum_e (b_e - (2 a_e + b_e)) = -2 sum_e a_e
    // assume 0 on v_1 and c on v2
    as := [-1/2 * e : e in edges];
    // sum_e -2a_e = v_2 = c
    c := -2 * &+as;
    bs := [c + e/2 : e in edges];
    cs := [0 : e in edges];
    assert forall{e : e in [1..#edges] | as[e] * 1^2 + bs[e] * 1 + cs[e] eq c};

    function values_with_subdivision(r)
        values := &join{{as[e] * (n/r)^2 + bs[e] * (n/r) + cs[e] : n in [0..r]} : e in [1..#edges]};
        return values;
    end function;

    printf "values with subdivision %o: %o\n", 1, values_with_subdivision(1);
    printf "values with subdivision %o: %o\n", 2, values_with_subdivision(2);

    printf "\n";
end for;