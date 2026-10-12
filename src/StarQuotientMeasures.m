intrinsic NTLRestrictMatrix(A::Mtrx, W::ModTupRng) -> Mtrx
{Return the matrix of the action of A restricted to the invariant subspace W.

INPUT:

- A -- square matrix acting on row vectors
- W -- subspace of the row space of A that is invariant under A

OUTPUT: the matrix of A on W with respect to Basis(W)
}
    B := Basis(W);
    TB := [b * A : b in B];
    return Matrix(BaseRing(A), [ Coordinates(W, v) : v in TB ]);
end intrinsic;

intrinsic NTLStarQuotientMeasures(p::RngIntElt, N::RngIntElt, r_list::SeqEnum) -> SeqEnum, SeqEnum, SeqEnum
{Return the Hecke matrices on the edges of the dual graph of X_0(pN)^* in characteristic p,
together with the stabiliser sizes and weights of the edges.

The edges of the dual graph of the semistable model of X_0(pN)^* in characteristic p
correspond to the supersingular points on X_0(pN) modulo the Atkin-Lehner involutions. The
i-th edge has length exactly stabilisers[i]*weights[i], so the semistable model of
X_0(pN)^* has stabilisers[i]*weights[i] - 1 extra components coming from the subdivision
of the i-th edge.

INPUT:

- p -- prime not dividing N
- N -- positive integer
- r_list -- sequence of integers coprime to pN

OUTPUT: a triple T, stabilisers, weights where

- T -- sequence of matrices; T[k] is the Hecke matrix T_r on the edges, for r = r_list[k]
- stabilisers -- sequence of integers; stabilisers[i] is the number of Atkin-Lehner
  involutions fixing the supersingular point on X_0(pN) corresponding to the i-th edge
- weights -- sequence of integers; weights[i] is half the number of automorphisms of the
  pair (E, G) representing the supersingular point on X_0(pN) corresponding to the i-th
  edge

WARNING:

This function is not implemented for all values of p, N and r_list; it is guaranteed to
raise an error if it hits an edge case it cannot handle.
}
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

    T := [NTLRestrictMatrix(HeckeOperator(B,r), Bstar) : r in r_list];
    return T, stabiliser_sizes, weights;
end intrinsic;

intrinsic NTLNrOfStarQuotientComponents(p::RngIntElt, N::RngIntElt) -> RngIntElt
{Return the number of components in characteristic p of the minimal regular model of
X_0(pN)^*.

INPUT:

- p -- prime not dividing N
- N -- positive integer

OUTPUT: a positive integer
}
    assert N mod p ne 0;
    _,stabiliser_sizes, weights := NTLStarQuotientMeasures(p, N, []);
    return 1 + &+([stabiliser_sizes[i]*weights[i]-1 : i in [1..#stabiliser_sizes]] cat [0]);
end intrinsic;

intrinsic NTLNrOfStarQuotientThickPoints(p::RngIntElt, N::RngIntElt) -> RngIntElt
{Return the number of supersingular points in characteristic p on X_0(pN)^* of thickness >
1.

This equals the number of loops of length > 1 in the dual graph in characteristic p of the
minimal regular model of X_0(pN)^*. This dual graph is one central vertex with a number of
loops attached to it.

INPUT:

- p -- prime not dividing N
- N -- positive integer

OUTPUT: a non-negative integer
}
    assert N mod p ne 0;
    _,stabiliser_sizes, weights := NTLStarQuotientMeasures(p, N, []);
    return #[i : i in [1..#stabiliser_sizes] | stabiliser_sizes[i]*weights[i] gt 1];
end intrinsic;

intrinsic NTLNrOfStarQuotientComponentConditions(N::RngIntElt) -> RngIntElt
{Return the total number of thick points of X_0(N)^* over all primes p dividing N.

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer

SEEALSO:

NTLNrOfStarQuotientThickPoints
}
    assert IsSquarefree(N);
    return &+[NTLNrOfStarQuotientThickPoints(p, N div p) : p in PrimeDivisors(N)];
end intrinsic;


intrinsic NTLStarQuotientGenus(N::RngIntElt) -> RngIntElt
{Return the genus of X_0(N)^*, computed from the dual graphs of its bad reductions.

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer

ALGORITHM:

Recursion on the number of prime divisors: for the smallest prime p | N the genus of
X_0(N)^* is the genus of X_0(N/p)^* plus the number of edges of the dual graph of X_0(N)^*
in characteristic p.

SEEALSO:

NTLGenusX0Nstar for a much faster closed formula.
}
    assert IsSquarefree(N);
    if N eq 1 then return 0; end if;
    p := PrimeDivisors(N)[1];
    g0 := NTLStarQuotientGenus(N div p);
    _,stabiliser_sizes, weights := NTLStarQuotientMeasures(p, N div p, []);
    return g0 + #stabiliser_sizes;
end intrinsic;

intrinsic NTLLocalHeightsAtp(p::RngIntElt, N::RngIntElt, rs::SeqEnum : warn := true) -> SeqEnum
{Return, for each r in rs, the set of values of the local height at p on the dual graph of
X_0(N)^* for the trace-zero part of T_r.

Intermediate results are printed.

INPUT:

- p -- prime exactly dividing N
- N -- squarefree positive integer
- rs -- sequence of integers coprime to N
- warn -- boolean (default: true); if true, raise an error explaining that this function
  is likely to be wrong

OUTPUT: a sequence of sets of rationals, one for each element of rs

WARNING:

This function is likely to have a bug: it computes traces of Hecke operators on the dual
graph of X_0(N)^* in characteristic p, but traces should be computed in characteristic 0,
as the dual graph in characteristic p only gives part of the trace.
}
    if warn then
        error "This function is likely to have a bug, as it computes traces of hecke operators on the dual graph of X_0(N)^* in characteristic p, but we should compute traces in characterstic 0, as the dual graph in characteristic p only give a part of the trace! Pass warn := false to ignore this warning.";
    end if;

    h_ps := [];
    Ts_raw, a, b := NTLStarQuotientMeasures(p, ExactQuotient(N, p), rs);
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
end intrinsic;




intrinsic NTLAtkinLehnerFixedSubspace(S::ModSym) -> ModTupFld
{Return the subspace of S fixed by all Atkin-Lehner operators.

INPUT:

- S -- space of weight-2 cuspidal modular symbols

OUTPUT: a subspace of VectorSpace(Rationals(), Dimension(S)), in the coordinates of S
}
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
end intrinsic;

intrinsic NTLStarQuotientCuspForms(N::RngIntElt) -> ModSym, ModTupFld
{Return the weight-2 cuspidal modular symbols of level N and sign +1, together with the
subspace fixed by all Atkin-Lehner operators.

The fixed subspace models S_2(Gamma_0(N))^W_N; in particular its dimension is the genus of
X_0(N)^*.

INPUT:

- N -- positive integer

OUTPUT: a pair S, Sstar where

- S -- the space of cuspidal modular symbols
- Sstar -- the Atkin-Lehner fixed subspace, as returned by NTLAtkinLehnerFixedSubspace
}
    S := CuspidalSubspace(ModularSymbols(N, 2, 1));
    return S, NTLAtkinLehnerFixedSubspace(S);
end intrinsic;

intrinsic NTLHeckeMatrixOnFixedSubspace(S::ModSym, Sstar::ModTupFld, r::RngIntElt) -> Mtrx
{Return the matrix of T_r on Sstar.

INPUT:

- S, Sstar -- as returned by NTLStarQuotientCuspForms
- r -- positive integer

OUTPUT: a square matrix over the rationals of size Dimension(Sstar)
}
    if Dimension(Sstar) eq 0 then
        return ZeroMatrix(Rationals(), 0, 0);
    end if;
    return NTLRestrictMatrix(HeckeOperator(S, r), Sstar);
end intrinsic;

intrinsic NTLStarQuotientHeckeMatrices(N::RngIntElt, r_list::SeqEnum) -> SeqEnum
{Return, for each r in r_list, the matrix of T_r acting on S_2(Gamma_0(N))^W_N.

INPUT:

- N -- positive integer
- r_list -- sequence of positive integers

OUTPUT: a sequence of square matrices of size g(X_0(N)^*)
}
    S, Sstar := NTLStarQuotientCuspForms(N);
    return [NTLHeckeMatrixOnFixedSubspace(S, Sstar, r) : r in r_list];
end intrinsic;

intrinsic NTLIsHeckeGenerator(T::Mtrx) -> BoolElt
{Return whether T generates the Hecke algebra.

Since the Hecke algebra tensor Q is a product of number fields, T_r generates it if and
only if its characteristic polynomial is squarefree (rem:Tp_generator in the paper).

INPUT:

- T -- the matrix of T_r on S_2(Gamma_0(N))^W_N

OUTPUT: a boolean
}
    return IsSquarefree(CharacteristicPolynomial(T));
end intrinsic;

intrinsic NTLStarQuotientHeightRelations(N::RngIntElt, r_list::SeqEnum) -> SeqEnum, SeqEnum, SeqEnum
{Return, for each r in r_list, the linear relations on a = (a_0, ..., a_(g-1)) for the
correspondence T_a = a_0 T_r^0 + ... + a_(g-1) T_r^(g-1) to have trivial height
contribution at all primes q | N.

These relations say that (T_a)_(j,j) = 0 for every j corresponding to a visible loop in
the dual graph of X_0(N)^* in characteristic q, for all q | N (Step 3 of
alg:clean_correspondences). All double points of X_0(N)^* in characteristic q are
non-split (lem:nonsplit), so a loop is visible if and only if its length is even
(def:visible_edge). Loops of odd length impose no condition
(lem:local_height_vanishes_long_edges).

INPUT:

- N -- squarefree positive integer
- r_list -- sequence of integers coprime to N

OUTPUT: a triple relations, is_hecke_generator, trace_relations where

- relations -- sequence of sequences of rows of length g; relations[i] contains the
  relations for T_r with r = r_list[i]. If it only contains zero rows, then any polynomial
  in T_r has trivial height contribution at all q | N
- is_hecke_generator -- sequence of booleans; whether T_r generates the Hecke algebra
- trace_relations -- sequence of sequences [Tr(T_r^0), ..., Tr(T_r^(g-1))], traces taken
  on S_2(Gamma_0(N))^W_N; T_a has trace 0 if and only if sum_k a_k Tr(T_r^k) = 0 (Step 1
  of alg:clean_correspondences)
}
    relations := [[] : i in [1..#r_list]];

    g := NTLStarQuotientGenus(N);

    // Whether T_r generates the Hecke algebra has to be decided on all of S_2(Gamma_0(N))^{W_N},
    // not on the dual graphs: otherwise polynomials f with f(T_r) = 0 can satisfy all relations.
    hecke_on_S := NTLStarQuotientHeckeMatrices(N, r_list);
    is_hecke_generator := [NTLIsHeckeGenerator(T) : T in hecke_on_S];
    assert forall{T : T in hecke_on_S | Nrows(T) eq g};
    trace_relations := [[Trace(T^k) : k in [0..g-1]] : T in hecke_on_S];

    for q in PrimeDivisors(N) do
        hecke, stabiliser_sizes, weights := NTLStarQuotientMeasures(q, N div q, r_list);

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
    return relations, is_hecke_generator, trace_relations;
end intrinsic;

intrinsic NTLStarQuotientDualGraph(N::RngIntElt)
{Print the lengths of the edges of the dual graph of X_0(N)^* in characteristic ell, for
every prime ell | N.

INPUT:

- N -- squarefree positive integer
}
    assert IsSquarefree(N);

    for ell in PrimeDivisors(N) do
        printf "Edges of the dual graph in characteristic ell = %o:\n", ell;
        _, stabiliser_sizes, weights := NTLStarQuotientMeasures(ell, N div ell, []);

        for j in [1..#stabiliser_sizes] do
            printf "Edge #%o of length %o.\n", j, stabiliser_sizes[j]*weights[j];
        end for;
    end for;
end intrinsic;

intrinsic NTLStarQuotientCleanCorrespondences(N::RngIntElt, r_list::SeqEnum) -> SeqEnum, SeqEnum, SeqEnum
{Return, for each r in r_list, a basis of the polynomials f in Z[x] such that the
correspondence T_f = f(T_r) has trivial height contribution at all primes q | N.

The T_f are not normalised to have trace 0. So a clean correspondence is only guaranteed
to exist if there are at least two polynomials for some r such that T_r generates the
Hecke algebra. If there is exactly one, there might still be a clean correspondence, but
one has to check directly that it has trace 0 on S_2(Gamma_0(N))^W_N.

The output can be passed to QCModAffine, since it takes suitable linear combinations to
ensure trace 0 automatically. At the time of writing QCModAffine raises an error if all
correspondences passed to it have trace 0.

INPUT:

- N -- squarefree positive integer
- r_list -- sequence of integers coprime to N

OUTPUT: a triple use_polys, is_hecke_generator, trace_relations where

- use_polys -- sequence of sequences of polynomials over the integers, one for each
  element of r_list
- is_hecke_generator, trace_relations -- as returned by NTLStarQuotientHeightRelations
}
    R<x> := PolynomialRing(Integers());
    use_polys := [[] : i in [1..#r_list]];
    relations, is_hecke_generator, trace_relations := NTLStarQuotientHeightRelations(N, r_list);
    for i in [1..#r_list] do
        if #relations[i] eq 0 then
            continue;
        end if;
        for v in Basis(Kernel(Transpose(Matrix(relations[i])))) do
            Append(~use_polys[i],&+[v[j]*x^(j-1) : j in [1..NumberOfColumns(v)]]);
        end for;
    end for;
    return use_polys, is_hecke_generator, trace_relations;
end intrinsic;

intrinsic NTLStarQuotientHasCleanCorrespondenceList(N::RngIntElt, r_list::SeqEnum) -> BoolElt, MonStgElt
{Return whether X_0(N)^* has at least one clean Hecke correspondence, using the first T_r
with r in r_list that generates the Hecke algebra.

When there are at least two kernel polynomials, a trace-zero linear combination always
exists, and as T_r generates the Hecke algebra, which has dimension g, it is non-zero.
When there is exactly one, it is checked directly whether it has trace 0 on
S_2(Gamma_0(N))^W_N.

INPUT:

- N -- squarefree positive integer
- r_list -- sequence of integers coprime to N

OUTPUT: a pair b, message where b is a boolean and message is a string explaining the
answer

An error is raised if no T_r with r in r_list generates the Hecke algebra.
}
    use_polys, is_hecke_generator, trace_relations := NTLStarQuotientCleanCorrespondences(N, r_list);
    hecke_gen_indices := [i : i in [1..#r_list] | is_hecke_generator[i]];
    error if #hecke_gen_indices eq 0,
        Sprintf("No T_r in r_list generates the Hecke algebra for N = %o.", N);
    i := hecke_gen_indices[1];
    if #use_polys[i] gt 1 then
        return true, Sprintf("There are at least 2 kernel polynomials for T_%o, so a trace-zero linear combination exists.", r_list[i]);
    elif #use_polys[i] eq 1 then
        r := r_list[i];
        f := use_polys[i][1];
        if &+[Rationals()| Coefficient(f, k)*trace_relations[i][k+1] : k in [0..#trace_relations[i]-1]] eq 0 then
            return true, Sprintf("There is exactly 1 kernel polynomial for T_%o, and it has trace 0 on the Atkin-Lehner +1 eigenspace of S_2(N).", r);
        else
            return false, Sprintf("There is exactly 1 kernel polynomial for T_%o, but it does not have trace 0 on the Atkin-Lehner +1 eigenspace of S_2(N).", r);
        end if;
    else
        return false, Sprintf("There are no kernel polynomials for T_%o, so there is no clean correspondence.", r_list[i]);
    end if;
end intrinsic;

intrinsic NTLStarQuotientHasCleanCorrespondence(N::RngIntElt) -> BoolElt, MonStgElt
{Return whether X_0(N)^* has at least one clean Hecke correspondence.

INPUT:

- N -- squarefree positive integer

OUTPUT: a pair b, message as returned by NTLStarQuotientHasCleanCorrespondenceList

ALGORITHM:

Find the smallest prime r coprime to N such that T_r generates the Hecke algebra (this
only needs modular symbols, which is cheap compared to the Brandt modules), and decide
using NTLStarQuotientHasCleanCorrespondenceList for that r.

An error is raised if no such r < 10000 exists.

TODO:

If g >= B + 2, where B is NTLNrOfStarQuotientComponentConditions, then the
(g-1)-dimensional space of non-scalar polynomials has dimension strictly greater than B,
so a clean correspondence exists by a dimension argument. This precheck is skipped for now
since NTLStarQuotientGenus is rather slow.
}
    // if NTLStarQuotientGenus(N) ge NTLNrOfStarQuotientComponentConditions(N) + 2 then
    //    return true;
    // end if;
    S, Sstar := NTLStarQuotientCuspForms(N);
    for r in [p : p in PrimesUpTo(10000) | not IsDivisibleBy(N, p)] do
        if NTLIsHeckeGenerator(NTLHeckeMatrixOnFixedSubspace(S, Sstar, r)) then
            return NTLStarQuotientHasCleanCorrespondenceList(N, [r]);
        end if;
    end for;
    error Sprintf("No Hecke generator found for N = %o among primes up to 10000.", N);
end intrinsic;


/*
p := 7;
N := p * 23;
r := 29;
h_p := NTLLocalHeightsAtp(p, N, [r]);
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
assert NTLStarQuotientHasCleanCorrespondence(246);
assert not NTLStarQuotientHasCleanCorrespondence(249);
*/
