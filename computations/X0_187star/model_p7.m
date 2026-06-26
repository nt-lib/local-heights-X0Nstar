SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";

X187star := X0Nstar(187);
pts := PointSearch(X187star, 30);

// Ambient projective space
P<X, Y, Z> := AmbientSpace(X187star);
 L1:= X + 2*Z;
 L2:= 2*X + Y - 2*Z;
 L3:= X + 2*Y + 2*Z;
 L4:= 2*X + Y - 2*Z;

S, Q, model_map, pts187, good_primes, inf_pts, bad_pts :=
                        find_and_test_model(L1, L2, L3, L4, X187star, x, y, pts
: printlevel:=1);

p := 7;
assert p in good_primes;

SetMemoryLimit(16 * 1024^3);

load "src/StarQuotientMeasures.m";
corrs := StarQuotientGoodCorrespondences(187, [p]);

polys := corrs[1];
good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p : use_polys := polys, number_of_correspondences := #polys-1, printlevel:=2 );




good_pts;
boo;
print "bad pts, fake pts, and bad disks:";
bad_pts;
fake_pts;
bad_disks;
