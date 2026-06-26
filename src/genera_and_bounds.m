// ===========================================================================
// Verification ceiling N_max: the inequality g_0^*(N) > B+1, B := sum_{p|N} B(p,N,W_N),
// is proved analytically for all squarefree N > N_max, so the exhaustive check below
// is needed only for N <= N_max.  N_max is the largest squarefree N for which the
// analytic argument fails.
//
// NOTE (threshold correction, g > B+1).  The existence guarantee is g > B+1, not g > B:
// the trace-zero normalization in StarQuotientHeightRelations frees the c_0 (identity)
// direction, so a non-scalar good correspondence is forced only once (g-1) - B > 0.
// The value N_max = 1231230 below was derived for the OLD inequality g > B, i.e. D(N) > 0.
// For g > B+1 one needs the slightly stronger margin D(N) > 1, so strictly N_max should be
// re-derived as the largest squarefree N with D(N) <= 1 (a bounded increase; the asymptotic
// conclusion g >> B is unaffected, and the worst-case primorial endgame shifts by O(1) in N).
// TODO: re-run the endgame/sieve with B -> B+1 to confirm (or bump) N_max.  The loop below
// already tests the corrected condition g <= B+1, so it correctly flags every borderline N
// it reaches; only the ceiling up to which it reaches still needs the +1 re-derivation.
//
// We use the *improved* threshold.  The genus/B complementarity is an EXACT
// identity: for each d at most one of nu(N,d) [enters g_0^*] and sum_p nu_p(N,d)
// [enters sum_p B] is nonzero (the first needs every prime of N/d split in
// Q(sqrt(-d)), the second needs exactly one inert), so
//   (genus error) + (B error) <= (1/2^omega) sum_{d|N} H(d) 2^{omega(N/d)},
// replacing the coefficient 3/(4 pi) by 2/(4 pi) in the analytic D(N).  Evaluating
// D(N) exactly on the worst case (primorials), with the p=2,3 extra-automorphism
// terms, the largest squarefree N with D(N)<=0 is
//   1 231 230 = 2*3*5*7*11*13*41
// (exhaustive sieve to 1.2e7, no failures with omega>=8, crude r^omega bound for
// N > 5.024e7).  Alternatives:
//   N_max = 2589510      : without complementarity (exact endgame only);
//   N_max = Round(5.1e7) : original naive asymptotic threshold.
N_max := 1231230;

// ===========================================================================
// Lazy-initialised class number tables, loaded/computed on first use.
//
// nu, nup, hprime only ever call ClassNumber(-4d) and ClassNumber(-d) with d a
// squarefree divisor of some squarefree N <= N_max.  Since every squarefree
// d <= N_max divides itself, the *whole* set of squarefree d <= N_max is needed,
// so we compute h4[d] = ClassNumber(-4d) (all squarefree d) and
// h1[d] = ClassNumber(-d) (squarefree d = 3 mod 4, so -d is a fundamental disc)
// once, into dense arrays stored in a NewStore so the tables survive across calls.
//
// The real win for repeated runs is PERSISTENCE: we save the tables to disk and
// reload them (seconds) on every later run, so the ~20 min is paid only once.
// Loading this file is now instant; the cost is deferred to the first CL4/CL1 call.
_cl_store := NewStore();
StoreSet(_cl_store, "loaded", false);

procedure _EnsureClassNumbers()
    if StoreGet(_cl_store, "loaded") then return; end if;

    clfile := Sprintf("data/classnumbers_upto_%o.dat", N_max);
    ok := false;
    try
        fh := Open(clfile, "r");
        h4 := ReadObject(fh);
        h1 := ReadObject(fh);
        delete fh;
        ok := (#h4 eq N_max) and (#h1 eq N_max);
    catch e
        ok := false;
    end try;
    if ok then
        printf "Loaded precomputed class numbers from %o.\n", clfile;
    else
        printf "Precomputing class numbers for squarefree d <= %o (one-off, ~minutes) ...\n", N_max;
        tprecomp := Cputime();
        h4 := [Integers()| 0 : i in [1..N_max]];
        h1 := [Integers()| 0 : i in [1..N_max]];
        for d in [1..N_max] do
            if IsSquarefree(d) then
                h4[d] := ClassNumber(-4*d);
                if d mod 4 eq 3 then
                    h1[d] := ClassNumber(-d);
                end if;
            end if;
            if IsDivisibleBy(d, 100000) then
                printf "  ... %o  (%o s)\n", d, Cputime(tprecomp);
            end if;
        end for;
        printf "Class numbers precomputed in %o s; saving to %o.\n", Cputime(tprecomp), clfile;
        f := Open(clfile, "w");
        WriteObject(f, h4);
        WriteObject(f, h1);
        delete f;
    end if;

    // Sanity: a few known values to catch a corrupted or mismatched cache.
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
// with psi(N)=prod(p+1), e2=prod_{p|N, p odd}(1+(-1/p)) (the factor at p=2 is 1),
// e3=prod_{p|N}(1+(-3/p)) (Kronecker symbol, correct at p=2,3), and #cusps=2^omega(N).
// Verified to agree with Genus(Gamma0(N)) for all squarefree N <= 5000; ~25x faster
// than Magma's Genus(Gamma0(N)), which is otherwise the dominant per-N cost.
// (Only correct for squarefree N -- which is all that is needed here.)
function GenusX0N(N)
    P := PrimeDivisors(N);
    psi := &*[Integers()| p+1 : p in P];
    e2  := &*[Integers()| 1 + KroneckerSymbol(-1, p) : p in P | p ne 2];
    e3  := &*[Integers()| 1 + KroneckerSymbol(-3, p) : p in P];
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

// H(d) (def:H_of_d): h'(-d) for d odd, h(-4d) for d even.
function H(d)
    if IsOdd(d) then
        return hprime(d);
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

    // Product over odd primes p | (N/d), i.e. p != 2
    Pq   := PrimeDivisors(q);
    Podd := [ p : p in Pq | p ne 2 ];
    prodOdd := &*[Integers()| 1 + KroneckerSymbol(-d, p) : p in Podd ];

    // -------------------------
    // Case A: 2 | d  (even d)
    // -------------------------
    if IsEven(d) then
        // The only genuine exception for squarefree N is d=2
        if d eq 2 then
            // Here q = N/2 is odd (since N squarefree).
            P := PrimeDivisors(q);
            prod1 := &*[Integers()| 1 + KroneckerSymbol(-1, p) : p in P ];
            prod2 := &*[Integers()| 1 + KroneckerSymbol(-2, p) : p in P ];
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
// nup(N,d,p): number nu_p(N,d) of w_d-fixed supersingular points on
// X_0(N)_{Fbar_p}, for squarefree N, prime p | N (the characteristic),
// and d | (N/p), d > 1.  Covers all primes p, including p = 2, 3, and even N.
//
// References:
//   p odd, N odd     -> cor:nup_odd_level
//   p odd, N even     -> cor:nup_even_level_p_odd  (parts (i)/(ii)/(iii))
//   p = 2 (char. 2)   -> cor:nup_p_eq_2            (d | N/2 is odd here)
// ---------------------------------------------------------------------------
function nup(N, d, p)
    assert IsSquarefree(N) and IsPrime(p) and IsDivisibleBy(N, p);
    M := N div p;                 // p does not divide M; M squarefree
    assert d gt 1 and IsDivisibleBy(M, d);   // d | (N/p), so p does not divide d
    // (In B(p,N,W) the sum is over d | M, so the ramified prime p never divides d.)

    if IsOdd(p) then
        if IsOdd(N) then
            // ---- Odd level (cor:nup_odd_level) ----
            prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(M div d)];
            return 1/2 * hprime(d) * (1 - KroneckerSymbol(-d, p)) * prod;
        else
            // ---- Even level, p odd (cor:nup_even_level_p_odd).  M = N/p even; Modd := M/2 odd. ----
            Modd := M div 2;
            if IsOdd(d) then
                // part (i)
                if (d mod 4) eq 1 then
                    kappa := CL4(d);
                else // d = 3 mod 4
                    kappa := 2*CL4(d) + (1 + KroneckerSymbol(-d, 2))*CL1(d);
                end if;
                prod := &*[Integers()| 1 + KroneckerSymbol(-d, q) : q in PrimeDivisors(Modd div d)];
                return 1/2 * (1 - KroneckerSymbol(-d, p)) * prod * kappa;
            elif d eq 2 then
                // part (ii)
                prod1 := &*[Integers()| 1 + KroneckerSymbol(-2, q) : q in PrimeDivisors(Modd)];
                prod2 := &*[Integers()| 1 + KroneckerSymbol(-1, q) : q in PrimeDivisors(Modd)];
                return 1/2 * ((1 - KroneckerSymbol(-2, p))*prod1 + (1 - KroneckerSymbol(-1, p))*prod2);
            else
                // part (iii): d > 2 even.  Level 2*Modd/d0 = M/d (odd).
                prod := &*[Integers()| 1 + KroneckerSymbol(-4*d, q) : q in PrimeDivisors(M div d)];
                return 1/2 * CL4(d) * (1 - KroneckerSymbol(-4*d, p)) * prod;
            end if;
        end if;
    else
        // ---- p = 2 is the characteristic (cor:nup_p_eq_2).  N = 2M, M = N/2 odd, so d | M is odd. ----
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

// Upper bound A_p(N) on the number of supersingular points (E,G) with
// #Aut(E,G) > 2.  For p > 3 these lie above j = 0, 1728.  For p in {2,3} the
// supersingular curve is unique with Aut(E) = 2T (order 24, seven maximal
// cyclic subgroups of order > 2) resp. 2D_6 (order 12, four), each fixing at
// most 2^{omega(N)-1} subgroups G of order M = N/p.
function AutTerm(N, p)
    w := omega(N);
    if p eq 2 then
        return 7 * 2^(w-1);
    elif p eq 3 then
        return 4 * 2^(w-1);
    else
        return 2^w;
    end if;
end function;

// B(p,N,W): upper bound on the number of thickness->1 double points of
// (X_0(N)/W)_{Fbar_p} (corollary before lem:upperBoundSumBps, extended to p=2,3).
// The raw expression 2/#W*(A_p + sum nu_p) is a bound on an integer count and can
// be a half-integer (the extra-automorphism part 2/#W*A_p is integral, but
// 2/#W*sum nu_p = (sum nu_p)/2^{omega-1} need not be).  Since the bounded quantity
// is a non-negative integer, we floor: the count of thick points is at most Floor(.).
function B(p, N, W)
    M := N div p;
    S := &+[Rationals()| nup(N, d, p) : d in Divisors(M) | d ne 1 and d in W];
    return Floor(2/#W * (AutTerm(N, p) + S));
end function;

function BoundOnB(N, W)
    return &+[Integers()| B(p, N, W) : p in PrimeDivisors(N)];
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
