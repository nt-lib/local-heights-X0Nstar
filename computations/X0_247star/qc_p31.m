// Quadratic Chabauty for X_0(247)^* at p = 31.
SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
AttachSpec("local_heights.spec");

N := 247;
X247star := X0Nstar(N);
_<W,X,Y,Z> := AmbientSpace(X247star);
pts := PointSearch(X247star, 30);

// The rational CM points found below are distinct and differ from the unique cusp of
// X_0(N)^*, so if they number one less than the rational points, all other points are CM.
AttachSpec("ShimuraCurveALQuotients/ShimuraQuotients.spec");
discs := RationalCMDiscs(CreateShimuraQuot(1, N, Set(Divisors(N))));
printf "Discriminants of the rational CM points, with the number of points: %o\n",
    [<d, discs[d]> : d in Reverse(Sort(Setseq(Keys(discs))))];
assert #pts eq 1 + &+[discs[d] : d in Keys(discs)];

// The plane model is the image of the map (x, y) = (L1/L2, L3/L4). The linear forms
// L1, ..., L4 were found by trying random linear forms with coefficients between -2 and 2
// until p was a good prime for the model.
L1 := -2*W;
L2 := -2*W - 2*X;
L3 := -Y - 2*Z;
L4 := 2*W + Z;
S, Q, model_map, model_pts, good_primes, inf_pts, bad_pts :=
    find_and_test_model(L1, L2, L3, L4, X247star, x, y, pts);
printf "Plane model Q(x, y) = 0 with\nQ = %o\n", Q;
printf "Good primes for the plane model: %o\n", good_primes;
printf "Known rational points on X_0(247)^* (PointSearch, height bound 30):\n%o\n", pts;
printf "Their images on the plane model:\n%o\n", model_pts;

p := 31;
assert p in good_primes;

use_polys := NTLStarQuotientCleanCorrespondences(N, [p])[1];
printf "Polynomials in T_p passed to QCModAffine (use_polys): %o\n", use_polys;
printf "Clean correspondences of trace 0, as polynomials in T_p: %o\n",
    NTLStarQuotientCleanCorrespondences(N, [p] : trace_zero := true)[1];

SetMemoryLimit(20 * 1024^3);
good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p : printlevel:=2,
    use_polys := use_polys, N:=20, number_of_correspondences := #use_polys - 1);

printf "Affine rational points on the plane model in good residue disks, found by quadratic Chabauty at p = %o:\n%o\n", p, good_pts;
assert boo;
assert #bad_pts eq 0;
assert #fake_pts eq 0;
assert #bad_disks eq 0;
// The points found by quadratic Chabauty are exactly the images of the known rational points.
assert forall{P : P in model_pts | P[3] eq 1};
assert Set(good_pts) eq {[P[1], P[2]] : P in model_pts};
printf "There are no bad points, fake points or bad disks, and the points found are exactly the images of the %o known rational points.\n", #pts;
