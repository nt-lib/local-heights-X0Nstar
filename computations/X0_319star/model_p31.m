SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";

X319star := X0Nstar(319);
_<W,X,Y,Z> := AmbientSpace(X319star);
pts := PointSearch(X319star, 30);
S, Q, model_map, pts319, good_primes, inf_pts, bad_pts := find_and_test_model(X, Z, W-X+Y, X, X319star, x, y, pts : max_prime:=100);
Q;
//y^3 + (-x^2 - 2*x + 12)*y^2 + (-2*x^5 - 5*x^4 - 8*x^3 - 15*x^2 - 10*x + 50)*y -
//    x^7 - 5*x^6 - 7*x^5 + x^4 - 11*x^3 - 47*x^2 - 5*x + 75

//  rec<recformat<p, infinite_points> |
//        p := 31,
//        infinite_points := [ (0 : 0 : 1 : 0), (29 : 29 : 1 : 0), (0 : 0 : 1 : 0)
//            ]>,
// no bad pts!

p := 31;

 polys := [x + 1, x^2 - 21, x^3 + 187];
good_pts, boo, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, p : printlevel:=2, use_polys := polys, N:=20, number_of_correspondences:=2);

good_pts;
//[
//    [ 0, -5 ],
//    [ 1, 0 ],
//    [ -1, -6 ]
//]
 boo;
 // true
bad_pts;
fake_pts;
bad_disks;
// all empty

// so we just need to find a model (with similar QC output)
// where the two F_31 points  (29 : 29 : 1 : 0), (0 : 0 : 1 : 0)
// are neither bad nor infinite