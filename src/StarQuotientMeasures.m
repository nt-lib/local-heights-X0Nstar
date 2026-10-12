function RestrictMatrix(A, W)
    B := Basis(W);
    TB := [b * A : b in B];
    return Matrix(BaseRing(A), [ Coordinates(W, v) : v in TB ]);
end function;

function StarQuotientMeasures(p, N, r_list)
    /* This function computes the Hecke matrix T_r on the edges of the dual graph of the semi stable model 
    of the star quotient X_0(pN)^* in characteristic p for every r in r_list. It returns a triple, T, stabilisers, weights . 
    
    - The list T contains the hecke matrices T_r for every r in r_list.
    - The stabilisers contain a list of how many Atkin-Lehner involutions fix the supersingular point on X_0(pN) corresponding to each edge
    - the weights list half of the number of automorphisms of the pair (E, G) representing the super singular point on X_0(pN)

    in particular the i-th edge should have length exactly stabilisers[i]*weights[i] . In particular the semi stable model of X_0(pN)^* will have stabilisers[i]*weights[i]-1 extra components coming from the subdivision of the i-th edge.

    Note that this function is not implemented for all values of p, N and r_list, it is guaranteed to raise an error if it hits an edge case it cannot handle.
    */
    assert GCD(p,N) eq 1 and forall{r : r in r_list | GCD(p*N,r) eq 1};
    B := BrandtModule(p,N);
    al_primes := AtkinLehnerPrimes(B);
    // need this since otherwise we cannot compute all the needed AL operators
    assert al_primes eq PrimeDivisors(Level(B));
    al_operators := [AtkinLehnerOperator(B,q) : q in al_primes];
    
    M := &* [w+1 : w in al_operators];
    Bstar := RowSpace(M);
    
    nonzero_entries := [[ i : i in [1..Ncols(v)] | v[i] ne 0 ] : v in Basis(Bstar)];
    stabiliser_sizes := [ExactQuotient(2^#al_primes, #entries) : entries in nonzero_entries];
    // I have a proof that the following assertion should be true in char >2, but can't hurt to
    // to double check. p=2, N=105 is an example where there is a stabiliser of size 4
    if p ne 2 then
        assert &and [size in [1,2] : size in stabiliser_sizes];
    end if;
    weights_ambient := MonodromyWeights(B);
    weights := [weights_ambient[entries[1]] : entries in nonzero_entries];

    T := [RestrictMatrix(HeckeOperator(B,r), Bstar) : r in r_list];
    return T, stabiliser_sizes, weights;
end function;

function NrOfStarQuotientComponents(p,N);
    // The number of components in characteristic p of the minimal regular model X_0(pN)^* 
    assert N mod p ne 0;
    _,stabiliser_sizes, weights := StarQuotientMeasures(p, N, []);
    return 1 + &+([stabiliser_sizes[i]*weights[i]-1 : i in [1..#stabiliser_sizes]] cat [0]);
end function;

function NrOfStarQuotientThickPoints(p,N);
    // The number of super singular points characteristic p on X_0(pN)^* of thickness > 1
    // This equals the number of loops of length in characteristic p in the dual graph of 
    // the minimal regular model X_0(N)^*
    // Note that the dual graph is one central vertex with a bunch of loops attached to it.
    assert N mod p ne 0;
    _,stabiliser_sizes, weights := StarQuotientMeasures(p, N, []);
    return #[i : i in [1..#stabiliser_sizes] | stabiliser_sizes[i]*weights[i] gt 1];
end function;

function NrOfStarQuotientComponentConditions(N);
    assert IsSquarefree(N);
    return &+[NrOfStarQuotientThickPoints(p, N div p) : p in PrimeDivisors(N)];
end function;


function StarQuotientGenus(N)
    assert IsSquarefree(N);
    if N eq 1 then return 0; end if;
    p := PrimeDivisors(N)[1];
    g0 := StarQuotientGenus(N div p);
    _,stabiliser_sizes, weights := StarQuotientMeasures(p, N div p, []);
    return g0 + #stabiliser_sizes;
end function;

function LocalHeightsAtp(p, N, rs : warn := true)
    if warn then
        error "This function is likely to have a bug, as it computes traces of hecke operators on the dual graph of X_0(N)^* in characteristic p, but we should compute traces in characterstic 0, as the dual graph in characteristic p only give a part of the trace! Pass warn := false to ignore this warning.";
    end if;

    h_ps := [];
    Ts_raw, a, b := StarQuotientMeasures(p, ExactQuotient(N, p), rs);
    // Trace-zero (primitive) part on the dual graph at p: subtract the scalar so that the
    // matrix has trace 0 on B* (dimension t_p = Nrows(Tr) = #loops = toric rank at p, NOT the
    // genus g and NOT the 2g-dim'l H^1_dR). With this factor the identity correspondence
    // (P_jj = 1, Trace = t_p) maps to 0, i.e. it has zero local height, as it must.
    Ts := [2*Nrows(Tr)*Tr - 2*Trace(Tr)*IdentityMatrix(Rationals(),Nrows(Tr)) : Tr in Ts_raw];
    traces := [Trace(Tr) : Tr in Ts_raw];
    printf "a = %o\nb = %o\n", a,b;
    for k -> r in rs do
        printf "r = %o\n", r;
        trace := traces[k];
        muT := [T[i,i] : i in [1..Nrows(T)]] where T := Ts[k];

        lengths := [a[i]*b[i] : i in [1..#a]];
        //mu_Ts := [[1/lengths[i] * T[i,i] : i in [1..Nrows(T)]] : T in Ts];
        //mu_Ts := [[T[i,i] : i in [1..Nrows(T)]] : T in Ts];
        printf "mu_T = %o\n", muT;

        // j_Gamma on edge[i] is a[i] * x^2 + b[i] * x
        // - 2 a [i] = muTs[i]
        // => a[i] = -1/2 mutTs[i]
        // a[i] * lengths[i]^2 + b[i] * lengths[i] = 0
        // => b[i] = - a[i] * lengths[i]

        as := [-1/2 * muT[j] : j in [1..#a]];
        bs := [-as[i] * lengths[i] : i in [1..#a]];
        assert forall{i : i in [1..#a] | 0 eq as[i] * lengths[i]^2 + bs[i] * lengths[i]};
        print &+bs - &+[2 * as[i] * lengths[i] + b[i] : i in [1..#a]], 1/2*trace;

        for i in [1..#a] do 
            //print [as[i] * (x/lengths[i])^2 + bs[i] * (x/lengths[i]) : x in [0..lengths[i]]];
            print [as[i] * x^2 + bs[i] * x : x in [0..lengths[i]]];
        end for;
        //h_p := &join{{as[i] * (x/lengths[i])^2 + bs[i] * (x/lengths[i]) : x in [0..lengths[i]]} 
        //                            : i in [1..#a]};
        h_p := &join{{as[i] * x^2 + bs[i] * x : x in [0..lengths[i]]} 
                                    : i in [1..#a]};
        Append(~h_ps, h_p);
    end for;
    return h_ps;
end function;




// Returns the subspace of S (weight-2 cuspidal modular symbols) fixed by all Atkin-Lehner
// operators, as a coordinate subspace of VectorSpace(Rationals(), Dimension(S)).
function AtkinLehnerFixedSubspace(S)
    N := Level(S);
    if Dimension(S) eq 0 then
        return VectorSpace(Rationals(), 0);
    end if;
    dim := Dimension(S);
    fixed := VectorSpace(Rationals(), dim);
    for pe in Factorization(N) do
        W := AtkinLehner(S, pe[1]^pe[2]);
        fixed := fixed meet Eigenspace(W, 1);
        if Dimension(fixed) eq 0 then break; end if;
    end for;
    return fixed;
end function;

function StarQuotientCuspForms(N)
    // Returns the weight-2 cuspidal modular symbols S of level N and sign +1, together with the
    // subspace Sstar fixed by all Atkin-Lehner operators, which models S_2(Gamma_0(N))^{W_N}.
    // In particular Dimension(Sstar) = g(X_0(N)^*).
    S := CuspidalSubspace(ModularSymbols(N, 2, 1));
    return S, AtkinLehnerFixedSubspace(S);
end function;

function HeckeMatrixOnFixedSubspace(S, Sstar, r)
    // The matrix of T_r on Sstar, for S, Sstar as returned by StarQuotientCuspForms.
    if Dimension(Sstar) eq 0 then
        return ZeroMatrix(Rationals(), 0, 0);
    end if;
    return RestrictMatrix(HeckeOperator(S, r), Sstar);
end function;

function StarQuotientHeckeMatrices(N, r_list)
    // Returns for each r in r_list the matrix of T_r acting on S_2(Gamma_0(N))^{W_N}.
    // These matrices have size g(X_0(N)^*).
    S, Sstar := StarQuotientCuspForms(N);
    return [HeckeMatrixOnFixedSubspace(S, Sstar, r) : r in r_list];
end function;

function IsHeckeGenerator(T)
    // T is the matrix of T_r on S_2(Gamma_0(N))^{W_N}. Since the Hecke algebra tensor Q is a
    // product of number fields, T_r generates it iff the characteristic polynomial of T is
    // squarefree (rem:Tp_generator in the paper).
    return IsSquarefree(CharacteristicPolynomial(T));
end function;

function StarQuotientHeightRelations(N, r_list)
    // Returns for each r in r_list the linear relations a :=(a_0,...,a_{g-1}) should satisfy in 
    // order for the correspondence given by:
    //    T_a := a_0(T_r)^0 + a_1(T_r)^1 + ... + a_{g-1}(T_r)^(g-1)
    // to have trivial height contribution at all primes q | N. The relations are given as a list 
    // of lists of length g, one for each r in r_list. The i-th list contains the relations for T_{r_i}.
    // If the list contains only zero rows, then there are no relations for T_{r_i} (i.e. any polynomial
    // in T_{r_i} has trivial height contribution at all q | N). The function also returns a list of
    // booleans indicating whether T_{r_i} generates the Hecke algebra (true) or not (false), and the
    // matrices of the T_{r_i} on S_2(Gamma_0(N))^{W_N}.
    //
    // These relations are just that T_{a,i,i} = 0 for all i corresponding to visible loops in the dual
    // graph of X_0(N)^* in characteristic q, for all q | N (Step 3 of alg:clean_correspondences).
    // All double points of X_0(N)^* in characteristic q are non-split (lem:nonsplit), so a loop is
    // visible iff its length is even (def:visible_edge). Loops of odd length impose no condition
    // (lem:local_height_vanishes_long_edges).
    relations := [[] : i in [1..#r_list]];

    g := StarQuotientGenus(N);

    // Whether T_r generates the Hecke algebra has to be decided on all of S_2(Gamma_0(N))^{W_N},
    // not on the dual graphs: otherwise polynomials f with f(T_r) = 0 can satisfy all relations.
    hecke_on_S := StarQuotientHeckeMatrices(N, r_list);
    is_hecke_generator := [IsHeckeGenerator(T) : T in hecke_on_S];
    assert forall{T : T in hecke_on_S | Nrows(T) eq g};

    for q in PrimeDivisors(N) do
        hecke, stabiliser_sizes, weights := StarQuotientMeasures(q, N div q, r_list);

        for i in [1..#r_list] do
            T := hecke[i]; 
            // Precompute the powers T^0,...,T^(g-1) once per (q, r); 
            powers := [T^k : k in [0..g-1]];
            for j in [1..#stabiliser_sizes] do
                // only visible loops, i.e. loops of even length, impose a condition
                if IsOdd(stabiliser_sizes[j]*weights[j]) then continue; end if;

                row := [powers[k+1][j,j] : k in [0..g-1]];
                Append(~relations[i], row);
            end for;
            if #relations[i] eq 0 then
                Append(~relations[i], [0 : k in [1..g]]);
            end if;
        end for;
    end for;
    return relations, is_hecke_generator, hecke_on_S;
end function;

procedure StarQuotientDualGraph(N)
    assert IsSquarefree(N);

    for ell in PrimeDivisors(N) do
        printf "Edges of the dual graph in characteristic ell = %o:\n", ell;
        _, stabiliser_sizes, weights := StarQuotientMeasures(ell, N div ell, []);

        for j in [1..#stabiliser_sizes] do
            printf "Edge #%o of length %o.\n", j, stabiliser_sizes[j]*weights[j];
        end for;
    end for;
end procedure;

function StarQuotientGoodCorrespondences(N, r_list)
    // Returns for each r in r_list a list of polynomials f(x) in Z[x] such that the correspondence
    // T_f := f(T_r) has trivial height contribution at all primes q | N. The function also returns 
    // a list of booleans indicating whether T_{r_i} generates the Hecke algebra (true) or not (false),
    // and the matrices of the T_{r_i} on S_2(Gamma_0(N))^{W_N}.
    // Note that these T_f have not yet been normalized to have trace 0. So an actualy good correspondence
    // is only guaranteed to exists if there are at least 2 polynomials in the list for some r
    // such that T_r generates the Hecke algebra.
    // If there is exactly 1, there might still be a good correspondence. But we need to check directly 
    // that it has trace 0 on the weight-2 cuspidal modular symbols of sign +1, restricted to the Atkin-Lehner +1 eigenspace.
    // 
    // Note that this is ok to pass to QCModAffine since QCModAffine will take suitable linear combinations to ensure trace 0 automatically.
    // Actually at the time of writing QCModAffine will raise an error if all correspondence you pass to it have trace 0.
    R<x> := PolynomialRing(Integers());
    use_polys := [[] : i in [1..#r_list]];
    relations, is_hecke_generator, hecke_on_S := StarQuotientHeightRelations(N, r_list);
    for i in [1..#r_list] do
        if #relations[i] eq 0 then
            continue;
        end if;
        for v in Basis(Kernel(Transpose(Matrix(relations[i])))) do
            Append(~use_polys[i],&+[v[j]*x^(j-1) : j in [1..NumberOfColumns(v)]]);
        end for;
    end for;
    return use_polys, is_hecke_generator, hecke_on_S;
end function;

function StarQuotientHasGoodCorrespondenceList(N, r_list)
    // Returns true iff X_0(N)^* has at least one good Hecke correspondence.
    // Raises an error if no T_r in r_list generates the Hecke algebra.
    // When there are >= 2 kernel polynomials a trace-zero linear combination always exists, and as
    // T_r generates the Hecke algebra, which has dimension g, it is non-zero.
    // When there is exactly 1, we check directly: f(T_r) has trace 0 on the weight-2
    // cuspidal modular symbols of sign +1, restricted to the Atkin-Lehner +1 eigenspace.
    use_polys, is_hecke_generator, hecke_on_S := StarQuotientGoodCorrespondences(N, r_list);
    hecke_gen_indices := [i : i in [1..#r_list] | is_hecke_generator[i]];
    error if #hecke_gen_indices eq 0,
        Sprintf("No T_r in r_list generates the Hecke algebra for N = %o.", N);
    i := hecke_gen_indices[1];
    if #use_polys[i] gt 1 then
        return true, Sprintf("There are at least 2 kernel polynomials for T_%o, so a trace-zero linear combination exists.", r_list[i]);
    elif #use_polys[i] eq 1 then
        r := r_list[i];
        f := use_polys[i][1];
        if Trace(Evaluate(f, hecke_on_S[i])) eq 0 then
            return true, Sprintf("There is exactly 1 kernel polynomial for T_%o, and it has trace 0 on the Atkin-Lehner +1 eigenspace of S_2(N).", r);
        else
            return false, Sprintf("There is exactly 1 kernel polynomial for T_%o, but it does not have trace 0 on the Atkin-Lehner +1 eigenspace of S_2(N).", r);
        end if;
    else
        return false, Sprintf("There are no kernel polynomials for T_%o, so there is no good correspondence.", r_list[i]);
    end if;
end function;

function StarQuotientHasGoodCorrespondence(N)
    // Returns true iff X_0(N)^* has at least one good Hecke correspondence.
    // Precheck: if g >= B + 2 (where B = NrOfStarQuotientComponentConditions) then the
    // (g-1)-dimensional space of non-scalar polynomials has dimension strictly greater than B,
    // so a good correspondence is guaranteed to exist by a dimension argument.
    // we skip it for now since StarQuotientGenus(N) is rather slow. Uncomment when we have a fast implementation of the genus.
    // if StarQuotientGenus(N) ge NrOfStarQuotientComponentConditions(N) + 2 then
    //    return true;
    // end if;
    // Otherwise find the smallest prime r coprime to N such that T_r generates the Hecke
    // algebra (this only needs modular symbols, which is cheap compared to the Brandt modules),
    // and decide using StarQuotientHasGoodCorrespondenceList for that r.
    S, Sstar := StarQuotientCuspForms(N);
    for r in [p : p in PrimesUpTo(10000) | not IsDivisibleBy(N, p)] do
        if IsHeckeGenerator(HeckeMatrixOnFixedSubspace(S, Sstar, r)) then
            return StarQuotientHasGoodCorrespondenceList(N, [r]);
        end if;
    end for;
    error Sprintf("No Hecke generator found for N = %o among primes up to 10000.", N);
end function;


/*
p := 7;
N := p * 23;
r := 29;
h_p := LocalHeightsAtp(p, N, [r]);
printf "local heights at %o: %o\n", p, h_p;

//N := 203;
//r := 19;
procedure values_of_local_heights(N, r)
    for p in PrimeDivisors(N) do
        h_p := local_heights_at_p(p, N, [r]);
        printf "local heights at %o: %o\n\n", p, h_p;
    end for;
end procedure;

values_of_local_heights(7*23, 29);

*/
/*
assert StarQuotientHasGoodCorrespondence(246);
assert not StarQuotientHasGoodCorrespondence(249);
*/
