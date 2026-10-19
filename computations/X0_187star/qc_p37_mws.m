SetMemoryLimit(40 * 1024^3);
AttachSpec("local_heights.spec");
SetPath("QCMod");
load "qc_modular.m";
load "divisor_heights.m";
load "qc_init_g2.m";
load "howe_zhu.m";
SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
load "src/misc.m";
load "src/mws_qc.m";

N := 187;
ALs := [11, 17];
seq_al := ALs;

X187star := X0NQuotient(N, ALs);
// Smooth, so no rational point below has several points above it on the desingularisation.
assert not IsSingular(X187star);
pts := PointSearch(X187star, 30);

// The rational CM points found below are distinct and differ from the unique cusp of
// X_0(N)^*, so if they number one less than the rational points, all other points are CM.
AttachSpec("ShimuraCurveALQuotients/ShimuraQuotients.spec");
discs := RationalCMDiscs(CreateShimuraQuot(1, N, Set(Divisors(N))));
printf "Discriminants of the rational CM points, with the number of points: %o\n",
    [<d, discs[d]> : d in Reverse(Sort(Setseq(Keys(discs))))];
assert #pts eq 1 + &+[discs[d] : d in Keys(discs)];

// Ambient projective space
P<X, Y, Z> := ProjectiveSpace(Rationals(), 2);

 L1:=X + Y + 2*Z;
 L2:=2*X - 2*Y + Z;
 L3:=Y - Z;
 L4:=-2*X + 2*Y - Z;

S, Q, model_map, pts187, good_primes, inf_pts, bad_pts :=
                        find_and_test_model(L1, L2, L3, L4, X187star, x, y, pts
: printlevel:=1);

primes := [37];
exponents := [Integers()| 1 : i in [1..#primes]];


J := JZero(N);
if #seq_al gt 0 then
    J_quo := Image(&*[1+AtkinLehnerOperator(J,i) : i in ALs]);
else
    J_quo := J;
end if;
J := J_quo;
gX := Dimension(J);
printf "genus = %o\n", gX;

ptsX := PointSearch(S, 30); // projective coordinates
index_of_base_point := [i : i in [1..#ptsX] | ptsX[i,3] ne 0][1]; // first index of ptsX such that the z-coordinate is non-zero
assert ptsX[index_of_base_point,3] ne 0;
base_pt := [ptsX[index_of_base_point,1]/ptsX[index_of_base_point,3], ptsX[index_of_base_point,2]/ptsX[index_of_base_point,3]]; 
printf "Using base point %o.\n", base_pt;

torsionMultiple := TorsionMultiple(J);
printf "#J(Q)_tors divides %o.\n", torsionMultiple;
assert torsionMultiple eq 1;

X := S;

pts := [ptsX[i] : i in [1..#ptsX] | i ne index_of_base_point];
base_point := ptsX[index_of_base_point];
// With D(x,y) := (x,y) - (-1/2,1/2), the known points satisfy D(-3,0) = 2D(-1/3,1) - 2D(2,1) and
// 2D(1/2,0) = -D(5/8,3/8) - 7D(2,1) + 6D(-1/3,1) (both divisors are principal), so D(2,1), D(-1/3,1),
// D(1/2,0) generate all D(P). The proof does not use this, only the checks d = 1 and filter_is_harmless.
// affine coordinates (x, y) = (X/Z, Y/Z) of a projective point
aff := func<P | [P[1]/P[3], P[2]/P[3]]>;
gens_xy := [[Rationals() | 2, 1], [Rationals() | -1/3, 1], [Rationals() | 1/2, 0]];
assert forall{a : a in gens_xy | exists{P : P in pts | aff(P) eq a}};
pts := [[P : P in pts | aff(P) eq a][1] : a in gens_xy] cat [P : P in pts | aff(P) notin gens_xy];
// coordinates of D(P) in the basis D(2,1), D(-1/3,1), D(1/2,0); checked against QC and in MWSieve
known_coords := [<[Rationals() | -1/2, 1/2], [0, 0, 0]>, <gens_xy[1], [1, 0, 0]>, <gens_xy[2], [0, 1, 0]>,
  <gens_xy[3], [0, 0, 1]>, <[Rationals() | -3, 0], [-2, 2, 0]>, <[Rationals() | 5/8, 3/8], [-7, 6, -2]>];
assert {aff(P) : P in ptsX} eq {t[1] : t in known_coords};
divisors := [* [* [[P[1]/P[3],P[2]/P[3]]], [[base_point[1]/base_point[3],base_point[2]/base_point[3]]] *]
                    : P in pts *];
assert gX eq 3;
splitting_indices := [ [1,0,0], [0,1,0], [0,0,1] ];

X := Curve(P, DefiningEquations(S)[1]);
printf "%o known points beside base point.\n", #pts;

exponents := [3];
qc_M := &*[primes[i]^exponents[i] : i in [1..#primes]];  // modulus
M := qc_M; // no 2-adic part: the d-check at 2 fails, so [J(Q) : G] might be even

// J(F_4127) is cyclic of order 67963257526 = 37^3 * 1341742

// primes p with 37 | #J(F_p) at which cosets are eliminated
mws_primes := [2, 29, 3617, 4127];

mws_primes := Sort(mws_primes);

bas := [[Eltseq(pt), Eltseq(base_point)] : pt in pts[1..gX]]; // same order as divisors[1..gX]
// 7 and 457 are bad primes of this plane model (X_0(187)^* itself has good reduction there)
sieve_primes := [v : v in mws_primes | v notin (PrimeDivisors(N) cat [7, 457])]; // the primes MWSieve uses
// the plane model has good reduction at the sieve primes
assert forall{v : v in sieve_primes | IsNonsingular(BaseChange(X, Bang(Rationals(), GF(v))))};
// d = 1 for rad(M) shows [J(Q) : G] is prime to M, hence G/MG = J(Q)/MJ(Q)
d := compute_d(X, bas, &*PrimeDivisors(M), sieve_primes);
printf "d = %o\n", d;
assert d eq 1;
assert filter_is_harmless(X, bas, sieve_primes, PrimeDivisors(M));

// quadratic Chabauty at p (Coleman data to precision N1), then the coefficients mod p^N2 of the fake
// and the known rational points in the basis divisors[1..gX]; retries with larger N1 on errors
function run_QC(p, N1, N2)
  assert N1 ge N2;

  try
    printf "Computing Coleman data with N1 = %o.\n", N1;
    data := coleman_data(Q, p, N1 : useU := false);

    printf "\nStarting quadratic Chabauty for p = %o.\n", p;
    R<x> := PolynomialRing(Rationals());
    corrs := NTLStarQuotientCleanCorrespondences(187, [p]);
    polys := corrs[1];
    printf "Polynomials in T_p passed to QCModAffine (use_polys): %o\n", polys;
    printf "Clean correspondences of trace 0, as polynomials in T_p: %o\n",
        NTLStarQuotientCleanCorrespondences(187, [p] : trace_zero := true)[1];
    time good_affine_rat_pts_xy, no_fake_pts, bad_affine_rat_pts_xy, data, fake_rat_pts, bad_Qppoints := QCModAffine(Q, p : use_polys := polys,
        number_of_correspondences := #polys-1, printlevel:=1, N := N1, prec := 40, base_point := base_pt);
    printf "There are %o fake rational points: %o\n", #fake_rat_pts, fake_rat_pts;
    printf "There are %o bad Q_%o-points.\n", #bad_Qppoints, p;

    // coefficients of the solutions in the basis divisors[1..gX], from Coleman data to precision N2;
    // retries with larger N2 (up to N1) on errors
    function recompute_Coleman_data(N1, N2, fake_rat_pts)
      try
        // Express the images of the solutions under Abel-Jacobi in terms of the generators mod p^N
        fake_rat_pts_new := [];
        for i in [1..#fake_rat_pts] do
          fake_rat_pts_new[i] := [ChangePrecision(fake_rat_pts[i,j], N2) : j in [1..2]];
          // lower precision for speed and to avoid issues in Coleman integrals.
        end for;
        printf "Re-computing Coleman data with N2 = %o.\n", N2;
        data := coleman_data(Q, p, N2 : useU := false);

        fake_coeffs_mod_pN, rat_coeffs_mod_pN := coefficients_mod_pN(fake_rat_pts_new, good_affine_rat_pts_xy, divisors, base_pt, splitting_indices, data : printlevel := 3); 
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

  return fake_coeffs_mod_pN, rat_coeffs_mod_pN, good_affine_rat_pts_xy;
end function;

fake_coeffs := [];
rat_coeffs := [];
for p in primes do
  fake_coeffs_mod_pN, rat_coeffs_mod_pN, rat_pts_xy := run_QC(p, 20, 8);
  // QC uses the basis bas: its coordinates of the known points agree with known_coords mod qc_M = 37^3
  assert {P : P in rat_pts_xy} eq {t[1] : t in known_coords};
  for i in [1..#rat_pts_xy] do
    coords := [t[2] : t in known_coords | t[1] eq rat_pts_xy[i]][1];
    assert forall{j : j in [1..gX] | IsDivisibleBy(rat_coeffs_mod_pN[i,j] - coords[j], qc_M)};
  end for;
  Append(~fake_coeffs, [fake_coeffs_mod_pN]);
  Append(~rat_coeffs, [rat_coeffs_mod_pN]);
end for;

print "done.";

printf "\nUse the Mordell-Weil sieve to show that the additional %o solutions aren't rational.\n", [#fake_coeffs_mod_pN[1] : fake_coeffs_mod_pN in fake_coeffs];
printf "generating cosets ...\n";
qc_fake_coeffs_mod_M := coeffs_CRT(fake_coeffs, primes, exponents);
printf "number of cosets: %o\n", #qc_fake_coeffs_mod_M;
qc_rat_coeffs_mod_M := [t[2] : t in known_coords];
printf "sieving with M = %o, S = %o\n", M, mws_primes;

printf "starting Mordell Weil sieve\n";
// SmoothBound exceeds every prime factor of #J(F_p), so MWSieve skips no prime
time done_fake := MWSieve(X, mws_primes, M, bas, base_point, qc_fake_coeffs_mod_M : 
  known_rat_coeffs := qc_rat_coeffs_mod_M, d := d, excluded := PrimeDivisors(N) cat [7, 457], SmoothBound := 77963257526, printlevel := 0);

printf "done with the MW sieve: %o\n", done_fake;
assert done_fake;
printf "No additional solutions are rational.\n";
