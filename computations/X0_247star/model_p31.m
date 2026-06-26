SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
X247star := X0Nstar(247);
P<W,X,Y,Z> := ProjectiveSpace(Rationals(), 3);
pts := PointSearch(X247star, 30);

L1 := -2*W;
L2 := -2*W - 2*X;
L3 := -Y - 2*Z;
L4 := 2*W + Z;
S, Q, model_map, pts201, good_primes, inf_pts, bad_pts :=
                        find_and_test_model(L1, L2, L3, L4, X247star, x, y,
pts);
good_primes;
Q;
p := good_primes[1];

load "src/StarQuotientMeasures.m";
corrs := StarQuotientGoodCorrespondences(247, [p]);

polys := corrs[1];


SetMemoryLimit(20 * 1024^3);



good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p : printlevel:=2, use_polys := polys, N:=20, number_of_correspondences:=2);

good_pts;
boo;
print "bad pts, fake pts, and bad disks:";
bad_pts;
fake_pts;
bad_disks;