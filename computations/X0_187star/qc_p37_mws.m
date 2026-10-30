SetMemoryLimit(40 * 1024^3);
AttachSpec("local_heights.spec");
SetPath("QCMod");
load "qc_modular.m";
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

p := 37;


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

qc_M := p^3;
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

fake_coeffs_mod_pN, rat_coeffs_mod_pN, rat_pts_xy, bad_Qppoints := run_QC(Q, N, p, 20, 8, divisors, base_pt, splitting_indices);
assert #bad_Qppoints eq 0;
// QC uses the basis bas: its coordinates of the known points agree with known_coords mod qc_M = 37^3
assert check_known_coords(rat_pts_xy, rat_coeffs_mod_pN, known_coords, qc_M);
print "done.";

printf "\nUse the Mordell-Weil sieve to show that the additional %o solutions aren't rational.\n", #fake_coeffs_mod_pN;
printf "sieving with M = %o, S = %o\n", M, mws_primes;

printf "starting Mordell Weil sieve\n";
time done_fake := MWSieve(X, mws_primes, M, bas, base_point, fake_coeffs_mod_pN :
  known_rat_coeffs := [t[2] : t in known_coords], d := d, excluded := PrimeDivisors(N) cat [7, 457]);

printf "done with the MW sieve: %o\n", done_fake;
assert done_fake;
printf "No additional solutions are rational.\n";
