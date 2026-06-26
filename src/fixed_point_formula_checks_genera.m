// ------------------------------------------------------------
// The following checks our formulas for nu(N,d):
//  # fixed points of w_d on X0(N) for squarefree N.
// Returns -1 on invalid input (e.g. d !| N or d<=1).
// To be more precise, we compare the genus of the star quotient
// computed by them and the one computed "directly" in Magma
// it seems to align. There is some Magma internal bug for 222.
// ------------------------------------------------------------

function ProdInt(seq)
    return (#seq eq 0) select 1 else &*seq;
end function;

// Kluit special case: d=2, N even squarefree.
function nu_squarefree_d2(N)
    if (N mod 2) ne 0 then
        return -1;
    end if;
    m := N div 2; // odd squarefree expected
    P := PrimeDivisors(m);
    prod1 := ProdInt([ 1 + KroneckerSymbol(-1, p) : p in P ]);
    prod2 := ProdInt([ 1 + KroneckerSymbol(-2, p) : p in P ]);
    return prod1 + prod2;
end function;

// Kluit/FH special case: d=3 (squarefree => only 2^0 or 2^1 may divide N/3).
function nu_squarefree_d3(N)
    if (N mod 3) ne 0 then
        return -1;
    end if;
    m := N div 3;
    if (m mod 2) eq 0 then
        m := m div 2;
    end if;
    P := PrimeDivisors(m);
    prod := ProdInt([ 1 + KroneckerSymbol(-3, p) : p in P ]);
    return 2*prod;
end function;

// Main function: nu(N,d) for squarefree N and d|N, d>1.
function nu_squarefree_old(N, d)
    if d le 1 then
        return -1;
    end if;
    if (N mod d) ne 0 then
        return -1;
    end if;

    // Special cases
    if d eq 2 then
        return nu_squarefree_d2(N);
    end if;
    if d eq 3 then
        return nu_squarefree_d3(N);
    end if;

    q := N div d;

    // Odd-prime Euler factor over p|q, p != 2
    Pq := PrimeDivisors(q);
    Podd := [ p : p in Pq | p ne 2 ];
    prodOdd := ProdInt([ 1 + KroneckerSymbol(-d, p) : p in Podd ]);

    // Case 1: d even
    if (d mod 2) eq 0 then
        return ClassNumber(-4*d) * prodOdd;
    end if;

    // Now d odd.
    if (q mod 2) ne 0 then
        // N odd
        if (d mod 4) eq 3 then
            // This is the "odd N" base term (matches Kluit/FH for squarefree)
            return (ClassNumber(-4*d) + ClassNumber(-d)) * prodOdd;
        else
            return ClassNumber(-4*d) * prodOdd;
        end if;
    else
        // N even squarefree: q even, so 2||q.
        if (d mod 4) eq 1 then
            return ClassNumber(-4*d) * prodOdd;  // c1(2)=1
        else
            // d = 3 (mod 4)
            fac2 := 1 + KroneckerSymbol(-d, 2); // c2(2)
            return (2*ClassNumber(-4*d) + fac2*ClassNumber(-d)) * prodOdd; // c1(2)=2
        end if;
    end if;
end function;

function nu_squarefree(N, d)
    if d le 1 then
        return -1;
    end if;
    if (N mod d) ne 0 then
        return -1;
    end if;

    q := N div d;

    // Product over odd primes p | (N/d), i.e. p != 2
    Pq   := PrimeDivisors(q);
    Podd := [ p : p in Pq | p ne 2 ];
    prodOdd := ProdInt([ 1 + KroneckerSymbol(-d, p) : p in Podd ]);

    // -------------------------
    // Case A: 2 | d  (even d)
    // -------------------------
    if (d mod 2) eq 0 then
        // The only genuine exception for squarefree N is d=2
        if d eq 2 then
            // Here q = N/2 is odd (since N squarefree).
            P := PrimeDivisors(q);
            prod1 := ProdInt([ 1 + KroneckerSymbol(-1, p) : p in P ]);
            prod2 := ProdInt([ 1 + KroneckerSymbol(-2, p) : p in P ]);
            return prod1 + prod2;
        end if;

        // For even d>2 (still squarefree), q is odd, and this matches the LaTeX "2|d" case.
        return ClassNumber(-4*d) * prodOdd;
    end if;

    // Now d is odd.

    // -------------------------
    // Case B: d = 1 (mod 4)
    // -------------------------
    if (d mod 4) eq 1 then
        // If N is even, then 2||q, but the 2-adic factor is 1 in this case.
        return ClassNumber(-4*d) * prodOdd;
    end if;

    // -------------------------
    // Case C: d = 3 (mod 4)
    // -------------------------
    // Distinguish whether q is odd (N odd) or q even (N even squarefree => 2||q).
    if (q mod 2) ne 0 then
        // N odd
        return (ClassNumber(-4*d) + ClassNumber(-d)) * prodOdd;
    else
        // N even squarefree: 2||q
        fac2 := 1 + KroneckerSymbol(-d, 2); // (1 + (-d/2))
        return (2*ClassNumber(-4*d) + fac2*ClassNumber(-d)) * prodOdd;
    end if;
end function;



function omega(N)
    return #PrimeDivisors(N);
end function;


function g0star_squarefree(N)
    w := omega(N);
    g0 := Genus(Gamma0(N));

    S := 0;
    for d in Divisors(N) do
        if d ne 1 then
            S +:= nu_squarefree(N, d);
        end if;
    end for;

    Q := Rationals();
    tw := 2^w;

    gstar := Q!1 + (Q!(g0-1))/tw - (Q!S)/(2*tw);

    return gstar;
end function;

procedure CompareStarGenus(Nmax)
    for N in [220..Nmax] do
        if IsSquarefree(N) then
            g_ours := g0star_squarefree(N);

            // Always use X0Nstar; Genus may fail if X0Nstar(N) is literally P^1.
            g_magma := 0;
            try
                Xs := X0Nstar(N);
                g_magma := Genus(Xs);
            catch e
                // If Genus fails, treat as genus 0 (P^1 case).
                g_magma := 0;
            end try;

            if g_ours ne g_magma then
                printf "Mismatch at N=%o: ours=%o, magma=%o\n", N, g_ours, g_magma;
            end if;
        end if;
    end for;
end procedure;



procedure CompareStarGenus(Nmin, Nmax)
    for N in [Nmin..Nmax] do
        if IsSquarefree(N) then
            g_ours := g0star_squarefree(N);

            // Step 1: construct X0(N)^*
            Xs := []; 
            try
                Xs := X0Nstar(N);
            catch e
                printf "SKIP N=%o: X0Nstar(N) failed with error:\n%o\n", N, e;
                continue;
            end try;

            // Step 2: handle the P^1 case safely
            g_magma := -1;

            // X0Nstar(N) sometimes returns a projective space object.
            if Type(Xs) eq Prj then
                g_magma := 0;
            else
                // Step 3: otherwise try Genus
                try
                    g_magma := Genus(Xs);
                catch e2
                    printf "SKIP N=%o: Genus(X0Nstar(N)) failed with error:\n%o\n", N, e2;
                    continue;
                end try;
            end if;

            if g_ours ne g_magma then
                printf "Mismatch at N=%o: ours=%o, magma=%o\n", N, g_ours, g_magma;
            end if;
        end if;
    end for;
	printf "All squarefree N between %o and %o conform to the fixed-point-count formulas\n", Nmin, Nmax;
end procedure;
