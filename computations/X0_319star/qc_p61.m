// Quadratic Chabauty for X_0(319)^* at p = 61.
SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
AttachSpec("local_heights.spec");

N := 319;
X319star := X0NQuotient(N, [11, 29]);
_<W,X,Y,Z> := AmbientSpace(X319star);
pts := PointSearch(X319star, 30);

// The plane model is the image of the map (x, y) = (L1/L2, L3/L4). The linear forms
// L1, ..., L4 were found by trying random linear forms with coefficients between -2 and 2
// until p was a good prime for the model.
L1 := W;
L2 := W + Y;
L3 := X;
L4 := Y + Z;
S, Q, model_map, model_pts, good_primes, inf_pts, bad_pts :=
    find_and_test_model(L1, L2, L3, L4, X319star, x, y, pts : max_prime := 100);
good_primes;
Q;

p := 61;
assert p in good_primes;

use_polys := NTLStarQuotientCleanCorrespondences(N, [p])[1];
printf "use_polys: %o\n", use_polys;

SetMemoryLimit(20 * 1024^3);
good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p : printlevel:=2,
    use_polys := use_polys, N:=20, number_of_correspondences := #use_polys - 1);

good_pts;
boo;
print "bad pts, fake pts, and bad disks:";
bad_pts;
fake_pts;
bad_disks;
assert boo;
assert #bad_pts eq 0;
assert #fake_pts eq 0;
assert #bad_disks eq 0;
