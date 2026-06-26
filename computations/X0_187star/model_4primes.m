SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
N := 187;
X187star := X0Nstar(N);
pts := PointSearch(X187star, 30);
// Ambient projective space
P<X, Y, Z> := ProjectiveSpace(Rationals(), 2);
L1 := X + Y;
L2 :=  X - Y - 2*Z;
L3 :=  -2*X - 2*Y - Z;
L4 := X - Y - 2*Z;
S, Q, model_map, pts187, good_primes, inf_pts, bad_pts :=
      find_and_test_model(L1, L2, L3, L4, X187star, x, y, pts : max_prime:=60, printlevel:=1);
Q;
A2<x,y> := AffineSpace(Rationals(), 2);
C := Curve(A2, y^4 + (44/5*x - 2)*y^3 + (573/20*x^2 - 25/2*x + 29/20)*y^2 +
(329/8*x^3 - 1033/40*x^2 + 231/40*x - 19/40)*y + 1761/80*x^4 - 177/10*x^3 + 229/40*x^2 - 9/10*x + 1/16);
print "Is this model singular:", IsSingular(C);
print "good primes for this model are ", good_primes;

p1 := good_primes[1];

load "src/StarQuotientMeasures.m";
corrs := StarQuotientGoodCorrespondences(N, [p1]);

polys := corrs[1];



good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p1 : printlevel:=2, use_polys := polys, N:=20, number_of_correspondences:=1);

good_pts;
boo;
print "bad pts, fake pts, and bad disks:";
bad_pts;
fake_pts;
bad_disks;
