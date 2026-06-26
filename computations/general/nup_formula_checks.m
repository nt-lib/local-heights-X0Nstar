// the code used to count nu_2(N, d)
// and check the conjectured formula
// (I've checked it for all odd squarefree M=N/d below 1000
// and a couple more


SetMemoryLimit(15 * 1024^3);
ZZ := Integers();

// modified class number h'(-d)
Hprime := function(d)
    assert d gt 0 and IsOdd(d);
    return (d mod 4 eq 1) select ClassNumber(-4*d) else (ClassNumber(-4*d) + ClassNumber(-d));
end function;

Prod := function(seq)
    return (#seq eq 0) select 1 else &*seq;
end function;

// (for) making sure we're not saved by rounding
ExactDiv := function(a, n)
    assert n gt 0;
    assert (a mod n) eq 0;
    return a div n;
end function;

// ------ predicted formula for nu_p(pM,d) ------

PredNu_pM := function(p, M, d)
    assert IsPrime(p);
    assert IsSquarefree(M) and (M mod p) ne 0;
    assert ((p*M) mod d) eq 0;

    // p = 2: N = 2M with M odd squarefree
    if p eq 2 then
        assert IsOdd(M);

        // d odd
        if IsOdd(d) then
            d0 := d;  assert (M mod d0) eq 0;
            primes := PrimeDivisors(M div d0);
            prod   := Prod([ 1 + KroneckerSymbol(-d0, q) : q in primes ]);

            if (d0 mod 4) eq 1 then
                return ExactDiv(ClassNumber(-4*d0) * prod, 2);
            end if;

            // d0 mod 4 = 3
            h4  := ClassNumber(-4*d0);
            hd  := ClassNumber(-d0);
            ks2 := KroneckerSymbol(-d0, 2); // ±1 for odd d0

            if d0 eq 3 then
                return ExactDiv((h4 + hd) * prod, 2);
            elif ks2 eq 1 then
                return 0;
            else
                return ExactDiv((h4 + hd) * prod, 4);
            end if;
        end if;

        // d even: d = 2*d0 with d0 | M
        d0 := d div 2;  assert d0 gt 0 and (M mod d0) eq 0;
        primes := PrimeDivisors(M div d0);

        // special: d = 2 (d0 = 1)
        if d0 eq 1 then
            prod8 := Prod([ 1 + KroneckerSymbol(-8, q) : q in primes ]);
            prod4 := Prod([ 1 + KroneckerSymbol(-4, q) : q in primes ]);
            return ExactDiv(ClassNumber(-8)*prod8 + ClassNumber(-4)*prod4, 2);
        end if;

        prod := Prod([ 1 + KroneckerSymbol(-8*d0, q) : q in primes ]);
        return ExactDiv(ClassNumber(-8*d0) * prod, 2);
    end if;

    // odd p branch (N = p*M is odd, hence d is odd)
    assert IsOdd(p*M) and IsOdd(d);

    d0 := (d mod p eq 0) select (d div p) else d;
    assert (M mod d0) eq 0;

    primes := PrimeDivisors(M div d0);
    prod   := Prod([ 1 + KroneckerSymbol(-d, q) : q in primes ]);

    num := Hprime(d) * (1 - KroneckerSymbol(-d, p)) * prod;
    return ExactDiv(num, 2);
end function;

// ----- brute fixed-point count on BrandtModule(p,M) ------

FixedPointsOfSignedPermMatrix := function(W)
    return #[ i : i in [1..Nrows(W)] | W[i,i] ne 0 ];
end function;

BruteNu_pM := function(p, M, d)
    assert IsPrime(p);
    assert IsSquarefree(M) and (M mod p) ne 0;
    assert ((p*M) mod d) eq 0;
    if p eq 2 then assert IsOdd(M); end if;

    BM := BrandtModule(p, M : ComputeGrams := true);
    n  := Dimension(BM);
    W  := IdentityMatrix(ZZ, n);

//    for q in PrimeDivisors(d) do
//        W *:= AtkinLehnerOperator(BM, q);
//    end for;

// documentation says that the 2nd argument has to be prime
// so we create W_d as product
// this works only for squarefree d 
    W := &*[ AtkinLehnerOperator(BM, q) : q in PrimeDivisors(d) ];

    return FixedPointsOfSignedPermMatrix(W);
end function;


Check_pM := function(p, M)
    assert IsPrime(p);
    assert IsSquarefree(M) and (M mod p) ne 0;
    N := p*M;

    print "--------------------------------------------------";
    print "Checking N = p*M with p =", p, " M =", M, " N =", N;

    divs := [ d : d in Divisors(N) | d gt 1 ];
    print "Divisors > 1:", divs;

    ok := true;
    for d in divs do
        brute := BruteNu_pM(p, M, d);
        pred  := PredNu_pM(p, M, d);
        print "d =", d, " brute =", brute, " pred =", pred;
        if brute ne pred then ok := false; end if;
    end for;

    print ok select "All d | N with d>1 match." else "Some mismatches occurred.";
    return ok;
end function;

CheckN := function(N)
    assert N gt 1 and IsSquarefree(N);

    ok := true;
    for p in PrimeDivisors(N) do
        M := N div p;

        if IsEven(M) then
            print "==================================================";
            print "Split: p =", p, " M =", M, " (so N =", N, ")";
            print "Skipped (out of scope: M even so N not squarefree).";
            continue;
        end if;

        print "==================================================";
        print "Split: p =", p, " M =", M, " (so N =", N, ")";

        if not Check_pM(p, M) then ok := false; end if;
    end for;

    print "==================================================";
    print ok select ("All tested splits passed for N = " cat IntegerToString(N))
             else ("Some tested split(s) failed for N = " cat IntegerToString(N));
    return ok;
end function;

print "ks2 =", KroneckerSymbol(-39, 2);
print "PredNu_pM(2,2145,39) =", PredNu_pM(2, 2145, 39);
print "BruteNu_pM(2,2145,39) =", BruteNu_pM(2, 2145, 39);

for M in [ n : n in [5..1000] | IsOdd(n) and IsSquarefree(n) ] do
    print "==== M =", M, "====";
    _ := Check_pM(2, M);
end for;