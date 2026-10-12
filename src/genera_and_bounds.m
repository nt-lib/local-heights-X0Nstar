/*
Auxiliary file for existence_of_relations.m: exact formulae for the fixed-point
counts nu(N,d) and nu_ell(N,d), the genus of X_0(N) and of its star quotient, and
the bound B(ell,N,W) on the number of loops of length > 1 in the dual graph at ell.
Class numbers h(-d), h(-4d) are read once from data/classnumbers_upto_<N_max>.dat.
N_max is the largest squarefree level with Phi(N) <= 0, so the analytic bound
ensures there's a fantastic correspondence for every N > N_max;
it is computed in final_inequality.m/txt.
*/

intrinsic NTLNMax() -> RngIntElt
{Return N_max, the largest squarefree level with Phi(N) <= 0.

For every N > N_max the analytic bound ensures there is a fantastic correspondence. It is
computed in computations/general/final_inequality.m. The class number tables cover all d
<= N_max.

OUTPUT: a positive integer
}
    return 2589510;
end intrinsic;

// ===========================================================================
// Lazy-initialised class number tables, loaded on first use.
//
// NTLnu, NTLnuell, NTLhprime only ever need h(-4d) and h(-d) with d a squarefree divisor of
// some squarefree N <= N_max.  Since every squarefree d <= N_max divides itself, the
// *whole* set of squarefree d <= N_max is needed, so we tabulate h4[d] = h(-4d)
// (all squarefree d) and h1[d] = h(-d) (squarefree d = 3 mod 4, so -d is a fundamental
// disc), and 0 otherwise, once, into dense arrays stored as attributes of Integers() so
// the tables survive across calls.
//
// The tables are stored in data/classnumbers_upto_<N_max>.dat, which is part of this
// repository; initialize_class_numbers.sh (optional) rewrites it from the LMFDB tables
// https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups (k = 0).
// Attaching the package is instant; the tables are read at the first NTLCL4/NTLCL1 call.
declare attributes RngInt: NTLClassNumbersH4, NTLClassNumbersH1;

procedure _EnsureClassNumbers()
    /*
    Read the class number tables from data/classnumbers_upto_<N_max>.dat and store
    them as attributes of Integers(), unless this has already been done.

    The file path is relative to the working directory, so Magma has to be run from the
    repository root. An error is raised if the file is missing or corrupted.
    */
    Z := Integers();
    if assigned Z`NTLClassNumbersH4 then return; end if;

    N_max := NTLNMax();
    clfile := Sprintf("data/classnumbers_upto_%o.dat", N_max);
    ok := true;
    try
        fh := Open(clfile, "r");
    catch e
        ok := false;
    end try;
    error if not ok, Sprintf("Cannot open %o; run initialize_class_numbers.sh to create it.", clfile);
    try
        h4 := ReadObject(fh);
        h1 := ReadObject(fh);
    catch e
        ok := false;
    end try;
    delete fh;
    ok := ok and Type(h4) eq SeqEnum and Type(h1) eq SeqEnum and #h4 eq N_max and #h1 eq N_max;
    error if not ok, Sprintf("%o was found but is not compatible with this version of Magma (or is corrupted); run initialize_class_numbers.sh to regenerate it.", clfile);
    printf "Loaded class numbers for d <= %o from %o.\n", N_max, clfile;

    // Sanity: a few known values to catch a corrupted or mismatched data file.
    assert h4[1] eq 1;    // h(-4)  = 1
    assert h4[5] eq 2;    // h(-20) = 2
    assert h1[3] eq 1;    // h(-3)  = 1
    assert h1[23] eq 3;   // h(-23) = 3

    Z`NTLClassNumbersH4 := h4;
    Z`NTLClassNumbersH1 := h1;
end procedure;

intrinsic NTLCL4(d::RngIntElt) -> RngIntElt
{Return the class number h(-4d) if d is squarefree, and 0 otherwise.

INPUT:

- d -- integer with 1 <= d <= N_max

OUTPUT: a non-negative integer

SEEALSO:

NTLNMax
}
    _EnsureClassNumbers();
    return Integers()`NTLClassNumbersH4[d];
end intrinsic;

intrinsic NTLCL1(d::RngIntElt) -> RngIntElt
{Return the class number h(-d) if d is squarefree and d = 3 mod 4, and 0 otherwise.

INPUT:

- d -- integer with 1 <= d <= N_max

OUTPUT: a non-negative integer

SEEALSO:

NTLNMax
}
    _EnsureClassNumbers();
    return Integers()`NTLClassNumbersH1[d];
end intrinsic;

// ===========================================================================

intrinsic NTLGenusX0N(N::RngIntElt) -> RngIntElt
{Return the genus of X_0(N) for squarefree N.

This is the closed form (rem:genusX0)

    g = 1 + psi(N)/12 - e_2/4 - e_3/3 - #cusps/2,

where

- psi(N) = prod_(ell | N) (ell + 1),
- e_2 = prod_(ell | N, ell odd) (1 + (-1/ell)); the factor at ell = 2 is 1,
- e_3 = prod_(ell | N) (1 + (-3/ell)), with the Kronecker symbol, which is also correct at
  ell = 2, 3,
- #cusps = 2^omega(N).

It has been verified to agree with Genus(Gamma0(N)) for all squarefree N <= 5000, and is
about 25 times faster than Magma's Genus(Gamma0(N)), which is otherwise the dominant per-N
cost.

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer

WARNING:

The result is only correct for squarefree N; this is not checked.
}
    P := PrimeDivisors(N);
    psi := &*[Integers()| ell+1 : ell in P];
    e2  := &*[Integers()| 1 + KroneckerSymbol(-1, ell) : ell in P | ell ne 2];
    e3  := &*[Integers()| 1 + KroneckerSymbol(-3, ell) : ell in P];
    return Integers()!(1 + psi/12 - e2/4 - e3/3 - 2^(#P)/2);
end intrinsic;

intrinsic NTLomega(N::RngIntElt) -> RngIntElt
{Return omega(N), the number of distinct prime divisors of N.

INPUT:

- N -- positive integer

OUTPUT: a non-negative integer
}
    return #PrimeDivisors(N);
end intrinsic;

intrinsic NTLhprime(d::RngIntElt) -> RngIntElt
{Return the modified class number h'(-d) (def:modified_classnumber).

This is h(-4d) if d = 1 mod 4 and h(-4d) + h(-d) if d = 3 mod 4.

INPUT:

- d -- odd squarefree positive integer

OUTPUT: a positive integer
}
    if d mod 4 eq 1 then
        return NTLCL4(d);
    elif d mod 4 eq 3 then
        return NTLCL1(d) + NTLCL4(d);
    end if;
end intrinsic;

intrinsic NTLH(d::RngIntElt) -> RngIntElt
{Return H(d) (def:H_of_d).

This is h'(-d) for odd d, h(-4d) for even d > 2, and H(2) = 2.

INPUT:

- d -- squarefree positive integer

OUTPUT: a positive integer

SEEALSO:

NTLhprime
}
    if IsOdd(d) then
        return NTLhprime(d);
    elif d eq 2 then
        return 2;   // h(-4) + h(-8), as in nu(N, 2)
    else
        return NTLCL4(d);
    end if;
end intrinsic;

intrinsic NTLnu(N::RngIntElt, d::RngIntElt) -> RngIntElt
{Return nu(N, d), the number of w_d-fixed points on X_0(N) in characteristic 0 (Kluit).

INPUT:

- N -- squarefree positive integer
- d -- divisor of N with d > 1

OUTPUT: a non-negative integer
}
    assert d ge 1 and IsDivisibleBy(N, d);

    q := N div d;

    // Product over odd primes ell | (N/d)
    Pq   := PrimeDivisors(q);
    Podd := [ ell : ell in Pq | ell ne 2 ];
    prodOdd := &*[Integers()| 1 + KroneckerSymbol(-d, ell) : ell in Podd ];

    // -------------------------
    // Case A: 2 | d  (even d)
    // -------------------------
    if IsEven(d) then
        // The only genuine exception for squarefree N is d=2
        if d eq 2 then
            // Here q = N/2 is odd (since N squarefree).
            P := PrimeDivisors(q);
            prod1 := &*[Integers()| 1 + KroneckerSymbol(-1, ell) : ell in P ];
            prod2 := &*[Integers()| 1 + KroneckerSymbol(-2, ell) : ell in P ];
            return prod1 + prod2;
        end if;

        // For even d>2 (still squarefree), q is odd, and this matches the LaTeX "2|d" case.
        return NTLCL4(d) * prodOdd;
    end if;

    // Now d is odd.

    // -------------------------
    // Case B: d = 1 (mod 4)
    // -------------------------
    if (d mod 4) eq 1 then
        // If N is even, then 2||q, but the 2-adic factor is 1 in this case.
        return NTLCL4(d) * prodOdd;
    end if;

    // -------------------------
    // Case C: d = 3 (mod 4)
    // -------------------------
    // Distinguish whether q is odd (N odd) or q even (N even squarefree => 2||q).
    if (q mod 2) ne 0 then
        // N odd
        return (NTLCL4(d) + NTLCL1(d)) * prodOdd;
    else
        // N even squarefree: 2||q
        fac2 := 1 + KroneckerSymbol(-d, 2); // (1 + (-d/2))
        return (2*NTLCL4(d) + fac2*NTLCL1(d)) * prodOdd;
    end if;
end intrinsic;

intrinsic NTLnuell(N::RngIntElt, d::RngIntElt, ell::RngIntElt) -> FldRatElt
{Return nu_ell(N, d), the number of w_d-fixed supersingular points on X_0(N) over Fbar_ell.

All primes ell are covered, including ell = 2, 3 and even N. The formulae are

- ell odd, N odd: cor:nuell_odd_level
- ell odd, N even: cor:nuell_even_level_ell_odd (parts (i), (ii), (iii))
- ell = 2: cor:nuell_ell_eq_2 (here d | N/2 is odd)

INPUT:

- N -- squarefree positive integer
- d -- divisor of N/ell with d > 1
- ell -- prime dividing N, the characteristic

OUTPUT: a non-negative integer, as an element of the rationals
}
    assert IsSquarefree(N) and IsPrime(ell) and IsDivisibleBy(N, ell);
    M := N div ell;                 // ell does not divide M; M squarefree
    assert d gt 1 and IsDivisibleBy(M, d);   // d | (N/ell), so ell does not divide d
    // (In B(ell,N,W) the sum is over d | M, so the ramified prime ell never divides d.)

    if IsOdd(ell) then
        if IsOdd(N) then
            // ---- Odd level (cor:nuell_odd_level) ----
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * NTLhprime(d) * (1 - KroneckerSymbol(-d, ell)) * prod;
        else
            // ---- Even level, ell odd (cor:nuell_even_level_ell_odd).  M = N/ell even; Modd := M/2 odd. ----
            Modd := M div 2;
            if IsOdd(d) then
                // part (i)
                if (d mod 4) eq 1 then
                    kappa := NTLCL4(d);
                else // d = 3 mod 4
                    kappa := 2*NTLCL4(d) + (1 + KroneckerSymbol(-d, 2))*NTLCL1(d);
                end if;
                prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(Modd div d)];
                return 1/2 * (1 - KroneckerSymbol(-d, ell)) * prod * kappa;
            elif d eq 2 then
                // part (ii)
                prod1 := &*[Integers()| 1 + KroneckerSymbol(-2, q) : q in PrimeDivisors(Modd)];
                prod2 := &*[Integers()| 1 + KroneckerSymbol(-1, q) : q in PrimeDivisors(Modd)];
                return 1/2 * ((1 - KroneckerSymbol(-2, ell))*prod1 + (1 - KroneckerSymbol(-1, ell))*prod2);
            else
                // part (iii): d > 2 even.  Level 2*Modd/d0 = M/d (odd).
                prod := &*[Integers()| 1 + KroneckerSymbol(-4*d, q) : q in PrimeDivisors(M div d)];
                return 1/2 * NTLCL4(d) * (1 - KroneckerSymbol(-4*d, ell)) * prod;
            end if;
        end if;
    else
        // ---- ell = 2 is the characteristic (cor:nuell_ell_eq_2).  N = 2M, M = N/2 odd, so d | M is odd. ----
        if (d mod 4) eq 1 then
            // part (i)
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * NTLCL4(d) * prod;
        else // d = 3 mod 4: part (ii); the factor at the ramified prime 2 is (1 - (-d/2))
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * NTLCL1(d) * (1 - KroneckerSymbol(-d, 2)) * prod;
        end if;
    end if;
end intrinsic;

intrinsic NTLAutTerm(N::RngIntElt, ell::RngIntElt) -> RngIntElt
{Return A_ell(N), an upper bound on the number of supersingular points (E, G) on X_0(N)
over Fbar_ell with #Aut(E, G) > 2.

For ell > 3 these lie above j = 0, 1728. For ell = 2, 3 the supersingular curve is unique
with Aut(E) = 2T (order 24, seven maximal cyclic subgroups of order > 2) resp. 2D_6 (order
12, four), each fixing at most 2^(omega(N)-1) subgroups G of order M = N/ell.

INPUT:

- N -- squarefree positive integer
- ell -- prime dividing N

OUTPUT: a positive integer
}
    w := NTLomega(N);
    if ell eq 2 then
        return 7 * 2^(w-1);
    elif ell eq 3 then
        return 4 * 2^(w-1);
    else
        return 2^w;
    end if;
end intrinsic;

intrinsic NTLB(ell::RngIntElt, N::RngIntElt, W::SeqEnum) -> RngIntElt
{Return B(ell, N, W), an upper bound on the number of double points of thickness > 1 of
X_0(N)/W over Fbar_ell (cor:Bpnw).

The raw expression 2/#W * (A_ell(N) + sum_d nu_ell(N, d)) bounds an integer count but can
be a half-integer, so the floor is returned.

INPUT:

- ell -- prime dividing N
- N -- squarefree positive integer
- W -- sequence of Hall divisors of N, representing a group of Atkin-Lehner involutions

OUTPUT: a non-negative integer

SEEALSO:

NTLAutTerm, NTLnuell
}
    M := N div ell;
    S := &+[Rationals()| NTLnuell(N, d, ell) : d in Divisors(M) | d ne 1 and d in W];
    return Floor(2/#W * (NTLAutTerm(N, ell) + S));
end intrinsic;

intrinsic NTLBoundOnB(N::RngIntElt, W::SeqEnum) -> RngIntElt
{Return sum_(ell | N) B(ell, N, W).

INPUT:

- N -- squarefree positive integer
- W -- sequence of Hall divisors of N, representing a group of Atkin-Lehner involutions

OUTPUT: a non-negative integer

SEEALSO:

NTLB
}
    return &+[Integers()| NTLB(ell, N, W) : ell in PrimeDivisors(N)];
end intrinsic;

intrinsic NTLGenusX0Nstar(N::RngIntElt) -> RngIntElt
{Return the genus of X_0(N)^* for squarefree N, using the Riemann-Hurwitz formula and the
fixed-point counts nu(N, d).

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer

SEEALSO:

NTLStarQuotientGenus for a (slow) independent computation.
}
    g := 1 + (NTLGenusX0N(N) - 1)/2^NTLomega(N) - 1/(2^(NTLomega(N)+1)) * &+[Integers()| NTLnu(N,d) : d in Divisors(N) | d ne 1];
    // The genus is an integer; coercing it asserts the rational formula came out integral
    // (a non-integer here would signal a bug in NTLGenusX0N or NTLnu).
    return Integers()!g;
end intrinsic;

intrinsic NTLHallDivisors(N::RngIntElt) -> SeqEnum
{Return the Hall divisors of N, i.e. the divisors d with gcd(d, N/d) = 1.

INPUT:

- N -- positive integer

OUTPUT: a sequence of positive integers
}
    return [d : d in Divisors(N) | Gcd(d, ExactQuotient(N,d)) eq 1];
end intrinsic;

// sanity test
/*
for N in [N : N in [1..350] | IsSquarefree(N)] do
    print N, NTLGenusX0Nstar(N), NTLStarQuotientGenus(N);
    assert NTLGenusX0Nstar(N) eq NTLStarQuotientGenus(N);
end for;*/
/*for N in [N : N in [50..100] | IsSquarefree(N)] do
    printf "N = %o:\n", N;
    WN := NTLHallDivisors(N);
    r_list := [p : p in PrimesUpTo(100) | not IsDivisibleBy(N, p)];
    for q in PrimeDivisors(N) do
        _, stabiliser_sizes, weights := NTLStarQuotientMeasures(q, N div q, r_list[1..1]);
        num_relations := 0;
        for j in [1..#stabiliser_sizes] do
            if stabiliser_sizes[j]*weights[j] eq 1 then continue; end if;
            num_relations +:= 1;
        end for;
        printf "relations for q = %o: %o <= %o\n", q, num_relations, NTLB(q,N,WN);
        assert num_relations le NTLB(q,N,WN);
    end for;
end for;*/
