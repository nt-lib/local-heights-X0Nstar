SetMemoryLimit(40 * 1024^3);
// Quadratic Chabauty at p = 47 and the Mordell-Weil sieve for X_0(430)^*. Run from the repo root.
AttachSpec("local_heights.spec");
SetPath("QCMod");
load "qc_modular.m";
SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
// loads QCMod again, so it must come before src/misc.m and src/mws_qc.m
load "ModularCurvesX0plusG4-6/preimage_checker.m";
load "src/misc.m";
load "src/mws_qc.m";

N := 430;
X430star := X0Nstar(N);
// the canonical model, a smooth plane quartic; the sieve works on it
assert not IsSingular(X430star);
pts := PointSearch(X430star, 30);
printf "Known rational points on X_0(430)^* (PointSearch, height bound 30): %o\n", pts;

AttachSpec("ShimuraCurveALQuotients/ShimuraQuotients.spec");
discs := RationalCMDiscs(CreateShimuraQuot(1, N, Set(Divisors(N))));
printf "Discriminants of the rational CM points, with the number of points: %o\n",
    [<d, discs[d]> : d in Reverse(Sort(Setseq(Keys(discs))))];
// So the 8 known points are the cusp, 6 rational CM points and one exceptional point, which was already found by
// Assaf, Hashimoto and Shnidman, "Exceptional points on Atkin-Lehner quotients", arXiv:2609.00516
// (https://github.com/sachihashimoto/AtkinLehnerQuotients).

// QC works with the plane model S below, which is singular but has the 8 known points in good affine residue
// discs at p = 47; its points at infinity lie over Q(sqrt(-2), i), which keeps QCModAffine's cup product cheap
P<X, Y, Z> := ProjectiveSpace(Rationals(), 2);
L1 := -2*X;
L2 := -2*Y;
L3 := 2*Z;
L4 := -2*X - 2*Y + Z;
S, Q, model_map, pts430, good_primes, inf_pts, bad_pts :=
    find_and_test_model(L1, L2, L3, L4, X430star, x, y, pts);
printf "Plane model Q(x, y) = 0 with\nQ = %o\n", Q;
// The plane model is only birational to X430star: there are no rational points where the map to the plane model
// is undefined, other than the known ones, and no unknown rational points mapping to the images of the known ones.
additional := additional_preimage(L1, L2, L3, L4, X430star, x, y, pts);
printf "additional preimages: %o\n", additional;
assert #additional eq 0;
p := 47;
assert p in good_primes;

J := JZero(N);
J := Image(&*[1 + AtkinLehnerOperator(J, q) : q in PrimeDivisors(N)]);
gX := Dimension(J);
printf "genus = %o\n", gX;
assert gX eq 3;
torsionMultiple := TorsionMultiple(J);
printf "#J(Q)_tors divides %o.\n", torsionMultiple;
assert torsionMultiple eq 1;

// affine coordinates (x, y) = (X/Z, Y/Z) of a projective point
aff := func<P | [P[1]/P[3], P[2]/P[3]]>;
// the known points on the model S; pts430[i] is the image of pts[i]
assert forall{T : T in pts430 | T[3] ne 0};
known_xy := [aff(T) : T in pts430];
// the point of X430star above a known point (x, y) of S
canon := func<T | pts[Index(known_xy, T)]>;
base_pt := [Rationals() | -1/2, 7/60];
// With D(x,y) := (x,y) - (-1/2,7/60), the other known points satisfy the relations in known_coords (all checked
// with IsPrincipal), so D(0,0), D(1,0), D(2/3,-1918/1215) generate all D(P). The proof only uses the checks below.
gens_xy := [[Rationals() | 0, 0], [Rationals() | 1, 0], [Rationals() | 2/3, -1918/1215]];
// coordinates of D(P) in the basis D(gens_xy); checked against QC and in MWSieve
known_coords := [<base_pt, [0, 0, 0]>, <gens_xy[1], [1, 0, 0]>, <gens_xy[2], [0, 1, 0]>,
  <gens_xy[3], [0, 0, 1]>, <[Rationals() | 0, 2/3], [-8, 7, 6]>, <[Rationals() | 1, 88/15], [-1, 2, 2]>,
  <[Rationals() | 2/3, 548/135], [-22, 17, 15]>, <[Rationals() | -1/2, 0], [-19, 15, 14]>];
assert {T : T in known_xy} eq {t[1] : t in known_coords};
divisors := [* [* [g], [base_pt] *] : g in gens_xy *]; // on S, for QC
splitting_indices := [ [1,0,0], [0,1,0], [0,0,1] ];
X := X430star; // the curve for the sieve
bas := [[Eltseq(canon(g)), Eltseq(canon(base_pt))] : g in gens_xy]; // same order as divisors
base_point := canon(base_pt);

exponents := [1];
qc_M := p^exponents[1];
// information mod 3 and 5 as well; not mod 2, since the d-check at 2 fails, so [J(Q) : G] might be even
aux_int := 3*5;
M := qc_M*aux_int;
// the good v < 400 with gcd(#J(F_v), M) > 1 at which cosets are eliminated, and 397, which eliminates none
// but is needed for d = 1 at 47
mws_primes := [13, 17, 23, 29, 59, 61, 89, 157, 163, 257, 263, 397];
sieve_primes := [v : v in mws_primes | v notin PrimeDivisors(N)]; // the primes MWSieve uses
// the canonical model has good reduction at the sieve primes
assert forall{v : v in sieve_primes | IsNonsingular(ChangeRing(X, GF(v)))};
// d = 1 for rad(M) shows [J(Q) : G] is prime to M, hence G/MG = J(Q)/MJ(Q)
d := compute_d(X, bas, &*PrimeDivisors(M), sieve_primes);
printf "d = %o\n", d;
assert d eq 1;
// MWSieve compares u-multiples, where u is the prime-to-M part of the LCM of the #J(F_v) (from the Hecke matrices
// on S_2(Gamma_0(N))^{W_N}), so only the parts for the primes of M matter; then the filter is harmless if
// [u J(F_v) : u red_v(G)] has no other prime factors
Ts := NTLStarQuotientHeckeMatrices(N, sieve_primes);
orders := [Integers()!Determinant((sieve_primes[i] + 1)*IdentityMatrix(Rationals(), gX) - Ts[i]) : i in [1..#sieve_primes]];
u := LCM([&*[Integers() | f[1]^f[2] : f in Factorization(n) | M mod f[1] ne 0] : n in orders]);
assert filter_is_harmless(X, bas, sieve_primes, PrimeDivisors(M) : mult := u);

// quadratic Chabauty at p (Coleman data to precision N1), then the coefficients mod p^N2 of the fake
// and the known rational points in the basis divisors[1..gX]; retries with larger N1 on errors
function run_QC(p, N1, N2)
  assert N1 ge N2;

  try
    printf "Computing Coleman data with N1 = %o.\n", N1;
    data := coleman_data(Q, p, N1 : useU := false);

    printf "\nStarting quadratic Chabauty for p = %o.\n", p;
    polys := NTLStarQuotientCleanCorrespondences(N, [p])[1];
    printf "Polynomials in T_p passed to QCModAffine (use_polys): %o\n", polys;
    printf "Clean correspondences of trace 0, as polynomials in T_p: %o\n",
        NTLStarQuotientCleanCorrespondences(N, [p] : trace_zero := true)[1];
    time good_affine_rat_pts_xy, no_fake_pts, bad_affine_rat_pts_xy, data, fake_rat_pts, bad_Qppoints := QCModAffine(Q, p : use_polys := polys,
        number_of_correspondences := #polys-1, printlevel:=1, N := N1, prec := 40, base_point := base_pt);
    printf "There are %o fake rational points.\n", #fake_rat_pts;
    printf "There are %o bad Q_%o-points.\n", #bad_Qppoints, p;

    // coefficients of the solutions in the basis divisors[1..gX], from Coleman data to precision N2;
    // retries with larger N2 (up to N1) on errors
    function recompute_Coleman_data(N1, N2, fake_rat_pts)
      try
        fake_rat_pts_new := [];
        for i in [1..#fake_rat_pts] do
          fake_rat_pts_new[i] := [ChangePrecision(fake_rat_pts[i,j], N2) : j in [1..2]];
        end for;
        printf "Re-computing Coleman data with N2 = %o.\n", N2;
        data := coleman_data(Q, p, N2 : useU := false);
        fake_coeffs_mod_pN, rat_coeffs_mod_pN := coefficients_mod_pN(fake_rat_pts_new, good_affine_rat_pts_xy, divisors, base_pt, splitting_indices, data);
      catch e
        printf "%o\n", e;
        if N2 lt N1 then
          N2 +:= 1;
        else
          error "reached N2 = N1";
        end if;
        printf "trying N1 = %o, N2 = %o ...\n", N1, N2;
        return $$(N1, N2, fake_rat_pts);
      end try;
      return fake_coeffs_mod_pN, rat_coeffs_mod_pN;
    end function;

    fake_coeffs_mod_pN, rat_coeffs_mod_pN := recompute_Coleman_data(N1, N2, fake_rat_pts);
  catch e
    N1 +:= 1;
    return $$(p, N1, N2);
  end try;

  return fake_coeffs_mod_pN, rat_coeffs_mod_pN, good_affine_rat_pts_xy, bad_Qppoints;
end function;

fake_coeffs_mod_pN, rat_coeffs_mod_pN, rat_pts_xy, bad_Qppoints := run_QC(p, 20, 8);
assert #bad_Qppoints eq 0;
// QC uses the basis bas: its coordinates of the known points agree with known_coords mod qc_M
assert {P : P in rat_pts_xy} eq {t[1] : t in known_coords};
for i in [1..#rat_pts_xy] do
  coords := [t[2] : t in known_coords | t[1] eq rat_pts_xy[i]][1];
  assert forall{j : j in [1..gX] | IsDivisibleBy(rat_coeffs_mod_pN[i,j] - coords[j], qc_M)};
end for;
print "done.";

printf "\nUse the Mordell-Weil sieve to show that the additional %o solutions aren't rational.\n", #fake_coeffs_mod_pN;
fake_coeffs_mod_M := combine_fake_good([[c mod qc_M : c in v] : v in fake_coeffs_mod_pN], qc_M, aux_int);
printf "adding information modulo %o: %o cosets\n", aux_int, #fake_coeffs_mod_M;
printf "sieving with M = %o, S = %o\n", M, sieve_primes;
time done_fake := MWSieve(X, mws_primes, M, bas, base_point, fake_coeffs_mod_M :
  known_rat_coeffs := [t[2] : t in known_coords], d := u, excluded := PrimeDivisors(N));

printf "done with the MW sieve: %o\n", done_fake;
assert done_fake;
printf "No additional solutions are rational.\n";
