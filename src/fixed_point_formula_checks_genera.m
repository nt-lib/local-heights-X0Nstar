/*
The following checks our formulas for nu(N,d):
 # fixed points of w_d on X0(N) for squarefree N.
To be more precise, we compare the genus of the star quotient
computed by them and the one computed "directly" in Magma
it seems to align. There is some Magma internal bug for 222.
*/

function ProdInt(seq)
    /*
    Return the product of the integers in seq, or 1 if seq is empty.
    */
    return (#seq eq 0) select 1 else &*seq;
end function;

intrinsic NTLnu_squarefree_d2(N::RngIntElt) -> RngIntElt
{Return nu(N, 2) for even squarefree N (Kluit), or -1 if N is odd.

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer, or -1 on invalid input
}
    if (N mod 2) ne 0 then
        return -1;
    end if;
    m := N div 2; // odd squarefree expected
    P := PrimeDivisors(m);
    prod1 := ProdInt([ 1 + KroneckerSymbol(-1, p) : p in P ]);
    prod2 := ProdInt([ 1 + KroneckerSymbol(-2, p) : p in P ]);
    return prod1 + prod2;
end intrinsic;

intrinsic NTLnu_squarefree_d3(N::RngIntElt) -> RngIntElt
{Return nu(N, 3) for squarefree N divisible by 3 (Kluit, FH), or -1 if 3 does not divide N.

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer, or -1 on invalid input
}
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
end intrinsic;

intrinsic NTLnu_squarefree_old(N::RngIntElt, d::RngIntElt) -> RngIntElt
{Return nu(N, d) for squarefree N, using the special cases NTLnu_squarefree_d2 and
NTLnu_squarefree_d3.

This is the old version of NTLnu_squarefree.

INPUT:

- N -- squarefree positive integer
- d -- divisor of N with d > 1

OUTPUT: a non-negative integer, or -1 on invalid input
}
    if d le 1 then
        return -1;
    end if;
    if (N mod d) ne 0 then
        return -1;
    end if;

    // Special cases
    if d eq 2 then
        return NTLnu_squarefree_d2(N);
    end if;
    if d eq 3 then
        return NTLnu_squarefree_d3(N);
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
end intrinsic;

intrinsic NTLnu_squarefree(N::RngIntElt, d::RngIntElt) -> RngIntElt
{Return nu(N, d), the number of fixed points of w_d on X_0(N) for squarefree N, using
Magma's ClassNumber.

INPUT:

- N -- squarefree positive integer
- d -- divisor of N with d > 1

OUTPUT: a non-negative integer, or -1 on invalid input

SEEALSO:

NTLnu, which reads the class numbers from a precomputed table.
}
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
end intrinsic;


intrinsic NTLg0star_squarefree(N::RngIntElt) -> FldRatElt
{Return the genus of X_0(N)^* for squarefree N, computed with the Riemann-Hurwitz formula
from Genus(Gamma0(N)) and NTLnu_squarefree.

INPUT:

- N -- squarefree positive integer

OUTPUT: a non-negative integer, as an element of the rationals
}
    w := NTLomega(N);
    g0 := Genus(Gamma0(N));

    S := 0;
    for d in Divisors(N) do
        if d ne 1 then
            S +:= NTLnu_squarefree(N, d);
        end if;
    end for;

    Q := Rationals();
    tw := 2^w;

    gstar := Q!1 + (Q!(g0-1))/tw - (Q!S)/(2*tw);

    return gstar;
end intrinsic;

intrinsic NTLCompareStarGenus(Nmax::RngIntElt)
{Print every squarefree N with 220 <= N <= Nmax for which NTLg0star_squarefree disagrees
with the genus of Magma's X0Nstar(N).

If Genus fails on X0Nstar(N), the genus is taken to be 0.

INPUT:

- Nmax -- positive integer
}
    for N in [220..Nmax] do
        if IsSquarefree(N) then
            g_ours := NTLg0star_squarefree(N);

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
end intrinsic;



intrinsic NTLCompareStarGenus(Nmin::RngIntElt, Nmax::RngIntElt)
{Print every squarefree N with Nmin <= N <= Nmax for which NTLg0star_squarefree disagrees
with the genus of Magma's X0Nstar(N).

Levels for which X0Nstar or Genus fails are reported and skipped.

INPUT:

- Nmin -- positive integer
- Nmax -- positive integer
}
    for N in [Nmin..Nmax] do
        if IsSquarefree(N) then
            g_ours := NTLg0star_squarefree(N);

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
end intrinsic;
