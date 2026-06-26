


SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";

X319star := X0NQuotient(319,[11, 29]);
_<W,X,Y,Z> := AmbientSpace(X319star);
X319star;
pts := PointSearch(X319star, 30);
S, Q, model_map, pts319, good_primes, inf_pts, bad_pts := find_and_test_model(W, W+Y, X, Y+Z, X319star, x, y, pts : max_prime := 100);
good_primes;
 // // [ 61 ]

p := 61;
polys := [x - 2, x^2 - 80, x^3 + 20];
good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p : printlevel:=2, use_polys := polys, N:=20, number_of_correspondences:=2);

good_pts;
boo;
print "bad pts, fake pts, and bad disks:";
bad_pts;
fake_pts;
bad_disks;