// Auxiliary file for existence_of_relations.m: exact formulae for the fixed-point
// counts nu(N,d) and nu_ell(N,d), the genus of X_0(N) and of its star quotient, and
// the bound B(ell,N,W) on the number of loops of length > 1 in the dual graph at ell.
// Class numbers h(-d), h(-4d) are read once from data/classnumbers_upto_<N_max>.dat.
// N_max is the largest squarefree level with Phi(N) <= 0, so the analytic bound
// ensures there's a fantastic correspondence for every N > N_max;
// it is computed in final_inequality.m/txt.

N_max := 2589510;

// ===========================================================================
// Lazy-initialised class number tables, loaded on first use.
//
// nu, nuell, hprime only ever need h(-4d) and h(-d) with d a squarefree divisor of
// some squarefree N <= N_max.  Since every squarefree d <= N_max divides itself, the
// *whole* set of squarefree d <= N_max is needed, so we tabulate h4[d] = h(-4d)
// (all squarefree d) and h1[d] = h(-d) (squarefree d = 3 mod 4, so -d is a fundamental
// disc), and 0 otherwise, once, into dense arrays stored in a NewStore so the tables
// survive across calls.
//
// The tables are stored in data/classnumbers_upto_<N_max>.dat, which is part of this
// repository; initialize_class_numbers.sh (optional) rewrites it from the LMFDB tables
// https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups (k = 0).
// Loading this file is instant; the tables are read at the first CL4/CL1 call.
_cl_store := NewStore();
StoreSet(_cl_store, "loaded", false);

procedure _EnsureClassNumbers()
    if StoreGet(_cl_store, "loaded") then return; end if;

    clfile := Sprintf("data/classnumbers_upto_%o.dat", N_max);
    ok := true;
    try
        fh := Open(clfile, "r");
    catch e
        ok := false;
    end try;
    error if not ok, Sprintf("Cannot open %o; run initialize_class_numbers.sh to create it.", clfile);
    h4 := ReadObject(fh);
    h1 := ReadObject(fh);
    delete fh;
    assert #h4 eq N_max and #h1 eq N_max;
    printf "Loaded class numbers for d <= %o from %o.\n", N_max, clfile;

    // Sanity: a few known values to catch a corrupted or mismatched data file.
    assert h4[1] eq 1;    // h(-4)  = 1
    assert h4[5] eq 2;    // h(-20) = 2
    assert h1[3] eq 1;    // h(-3)  = 1
    assert h1[23] eq 3;   // h(-23) = 3

    StoreSet(_cl_store, "h4", h4);
    StoreSet(_cl_store, "h1", h1);
    StoreSet(_cl_store, "loaded", true);
end procedure;

function CL4(d)
    _EnsureClassNumbers();
    return StoreGet(_cl_store, "h4")[d];
end function;

function CL1(d)
    _EnsureClassNumbers();
    return StoreGet(_cl_store, "h1")[d];
end function;

// ===========================================================================

// Closed form for the genus of X_0(N), N squarefree (rem:genusX0):
//   g = 1 + psi(N)/12 - e2/4 - e3/3 - (#cusps)/2,
// with psi(N)=prod(ell+1), e2=prod_{ell|N, ell odd}(1+(-1/ell)) (the factor at ell=2 is 1),
// e3=prod_{ell|N}(1+(-3/ell)) (Kronecker symbol, correct at ell=2,3), and #cusps=2^omega(N).
// Verified to agree with Genus(Gamma0(N)) for all squarefree N <= 5000; ~25x faster
// than Magma's Genus(Gamma0(N)), which is otherwise the dominant per-N cost.
// (Only correct for squarefree N -- which is all that is needed here.)
function GenusX0N(N)
    P := PrimeDivisors(N);
    psi := &*[Integers()| ell+1 : ell in P];
    e2  := &*[Integers()| 1 + KroneckerSymbol(-1, ell) : ell in P | ell ne 2];
    e3  := &*[Integers()| 1 + KroneckerSymbol(-3, ell) : ell in P];
    return Integers()!(1 + psi/12 - e2/4 - e3/3 - 2^(#P)/2);
end function;

function omega(N)
    return #PrimeDivisors(N);
end function;

// hprime(d) = h'(-d) for d>0 odd squarefree (def:modified_classnumber):
//   h(-4d)            if d = 1 mod 4,
//   h(-4d) + h(-d)    if d = 3 mod 4.
function hprime(d)
    if d mod 4 eq 1 then
        return CL4(d);
    elif d mod 4 eq 3 then
        return CL1(d) + CL4(d);
    end if;
end function;

// H(d) (def:H_of_d): h'(-d) for d odd, h(-4d) for even d > 2, and H(2) = 2.
function H(d)
    if IsOdd(d) then
        return hprime(d);
    elif d eq 2 then
        return 2;   // h(-4) + h(-8), as in nu(N, 2)
    else
        return CL4(d);
    end if;
end function;

// ---------------------------------------------------------------------------
// nu(N,d): number of w_d-fixed points on X_0(N) in characteristic 0 (Kluit),
// for squarefree N and d | N, d > 1.  Already valid for all squarefree N.
// ---------------------------------------------------------------------------
function nu(N, d)
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
        return CL4(d) * prodOdd;
    end if;

    // Now d is odd.

    // -------------------------
    // Case B: d = 1 (mod 4)
    // -------------------------
    if (d mod 4) eq 1 then
        // If N is even, then 2||q, but the 2-adic factor is 1 in this case.
        return CL4(d) * prodOdd;
    end if;

    // -------------------------
    // Case C: d = 3 (mod 4)
    // -------------------------
    // Distinguish whether q is odd (N odd) or q even (N even squarefree => 2||q).
    if (q mod 2) ne 0 then
        // N odd
        return (CL4(d) + CL1(d)) * prodOdd;
    else
        // N even squarefree: 2||q
        fac2 := 1 + KroneckerSymbol(-d, 2); // (1 + (-d/2))
        return (2*CL4(d) + fac2*CL1(d)) * prodOdd;
    end if;
end function;

// ---------------------------------------------------------------------------
// nuell(N,d,ell): number nu_ell(N,d) of w_d-fixed supersingular points on
// X_0(N)_{Fbar_ell}, for squarefree N, prime ell | N (the characteristic),
// and d | (N/ell), d > 1.  Covers all primes ell, including ell = 2, 3, and even N.
//
// References:
//   ell odd, N odd     -> cor:nuell_odd_level
//   ell odd, N even     -> cor:nuell_even_level_ell_odd  (parts (i)/(ii)/(iii))
//   ell = 2 (char. 2)   -> cor:nuell_ell_eq_2            (d | N/2 is odd here)
// ---------------------------------------------------------------------------
function nuell(N, d, ell)
    assert IsSquarefree(N) and IsPrime(ell) and IsDivisibleBy(N, ell);
    M := N div ell;                 // ell does not divide M; M squarefree
    assert d gt 1 and IsDivisibleBy(M, d);   // d | (N/ell), so ell does not divide d
    // (In B(ell,N,W) the sum is over d | M, so the ramified prime ell never divides d.)

    if IsOdd(ell) then
        if IsOdd(N) then
            // ---- Odd level (cor:nuell_odd_level) ----
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * hprime(d) * (1 - KroneckerSymbol(-d, ell)) * prod;
        else
            // ---- Even level, ell odd (cor:nuell_even_level_ell_odd).  M = N/ell even; Modd := M/2 odd. ----
            Modd := M div 2;
            if IsOdd(d) then
                // part (i)
                if (d mod 4) eq 1 then
                    kappa := CL4(d);
                else // d = 3 mod 4
                    kappa := 2*CL4(d) + (1 + KroneckerSymbol(-d, 2))*CL1(d);
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
                return 1/2 * CL4(d) * (1 - KroneckerSymbol(-4*d, ell)) * prod;
            end if;
        end if;
    else
        // ---- ell = 2 is the characteristic (cor:nuell_ell_eq_2).  N = 2M, M = N/2 odd, so d | M is odd. ----
        if (d mod 4) eq 1 then
            // part (i)
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * CL4(d) * prod;
        else // d = 3 mod 4: part (ii); the factor at the ramified prime 2 is (1 - (-d/2))
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * CL1(d) * (1 - KroneckerSymbol(-d, 2)) * prod;
        end if;
    end if;
end function;

// Upper bound A_ell(N) on the number of supersingular points (E,G) with
// #Aut(E,G) > 2.  For ell > 3 these lie above j = 0, 1728.  For ell in {2,3} the
// supersingular curve is unique with Aut(E) = 2T (order 24, seven maximal
// cyclic subgroups of order > 2) resp. 2D_6 (order 12, four), each fixing at
// most 2^{omega(N)-1} subgroups G of order M = N/ell.
function AutTerm(N, ell)
    w := omega(N);
    if ell eq 2 then
        return 7 * 2^(w-1);
    elif ell eq 3 then
        return 4 * 2^(w-1);
    else
        return 2^w;
    end if;
end function;

// B(ell,N,W): upper bound on the number of thickness > 1 double points of
// (X_0(N)/W)_{Fbar_ell} (corollary cor:Bpnw).
// The raw expression 2/#W*(A_ell + sum nu_ell) is a bound on an integer count and can
// be a half-integer so the count of thick points is at most Floor(.).
function B(ell, N, W)
    M := N div ell;
    S := &+[Rationals()| nuell(N, d, ell) : d in Divisors(M) | d ne 1 and d in W];
    return Floor(2/#W * (AutTerm(N, ell) + S));
end function;

function BoundOnB(N, W)
    return &+[Integers()| B(ell, N, W) : ell in PrimeDivisors(N)];
end function;

function GenusX0Nstar(N)
    g := 1 + (GenusX0N(N) - 1)/2^omega(N) - 1/(2^(omega(N)+1)) * &+[Integers()| nu(N,d) : d in Divisors(N) | d ne 1];
    // The genus is an integer; coercing it asserts the rational formula came out integral
    // (a non-integer here would signal a bug in GenusX0N or nu).
    return Integers()!g;
end function;

function HallDivisors(N)
    return [d : d in Divisors(N) | Gcd(d, ExactQuotient(N,d)) eq 1];
end function;

// sanity test
/*
for N in [N : N in [1..350] | IsSquarefree(N)] do
    print N, GenusX0Nstar(N), StarQuotientGenus(N);
    assert GenusX0Nstar(N) eq StarQuotientGenus(N);
end for;*/
/*for N in [N : N in [50..100] | IsSquarefree(N)] do
    printf "N = %o:\n", N;
    WN := HallDivisors(N);
    r_list := [p : p in PrimesUpTo(100) | not IsDivisibleBy(N, p)];
    for q in PrimeDivisors(N) do
        _, stabiliser_sizes, weights := StarQuotientMeasures(q, N div q, r_list[1..1]);
        num_relations := 0;
        for j in [1..#stabiliser_sizes] do
            if stabiliser_sizes[j]*weights[j] eq 1 then continue; end if;
            num_relations +:= 1;
        end for;
        printf "relations for q = %o: %o <= %o\n", q, num_relations, B(q,N,WN);
        assert num_relations le B(q,N,WN);
    end for;
end for;*/
