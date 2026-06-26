SetDebugOnError(true);
/*X := X0Nstar(319);
pts := PointSearch(X,100);
for pt in pts do
    X2 := ProjectionFromNonsingularPoint(X,pt);
    P2<x,y,z> := AmbientSpace(X2);
    print X2;
end for; */

_<x> := PolynomialRing(Rationals());

use_polys := [
<2, [x, x^2 - 2, x^3 + 2]>,
<3, [x + 1, x^2 - 3, x^3 + 5]>,
<5, [x + 1, x^2 - 3, x^3 + 9]>,
<7, [x, x^2 - 4, x^3]>,
<13, [x + 2, x^2 - 22, x^3 + 120]>,
<17, [x, x^2 - 26, x^3 + 52]>,
<19, [x, x^2 - 2, x^3 + 2]>,
<23, [x + 1, x^2 - 13, x^3 + 33]>,
<31, [x + 1, x^2 - 21, x^3 + 187]>,
<37, [x + 1, x^2 - 19, x^3 + 109]>,
<41, [x + 10, x^2 - 112, x^3 + 1288]>,
<43, [x - 4, x^2 - 86, x^3 - 666]>,
<47, [x + 2, x^2 - 22, x^3 + 190]>,
<53, [x - 6, x^2 - 86, x^3 - 1060]>,
<59, [x + 5, x^2 - 51, x^3 + 553]>,
<61, [x - 2, x^2 - 80, x^3 + 20]>,
<67, [x + 1, x^2 - 77, x^3 - 49]>,
<71, [x + 3, x^2 - 79, x^3 + 743]>,
<73, [x + 4, x^2 - 90, x^3 + 518]>,
<79, [x + 2, x^2 - 58, x^3 + 386]>,
<83, [x + 2, x^2 - 134, x^3 + 156]>,
<89, [x + 7, x^2 - 135, x^3 + 1287]>,
<97, [x - 1, x^2 - 107, x^3 - 485]>,
<101, [x + 6, x^2 - 90, x^3 + 702]>,
<103, [x - 4, x^2 - 116, x^3 - 1520]>,
<107, [x - 4, x^2 - 74, x^3 - 710]>,
<109, [x, x^2 - 182, x^3 - 650]>,
<113, [x + 7, x^2 - 85, x^3 + 1315]>,
<127, [x - 6, x^2 - 110, x^3 - 1114]>,
<131, [x + 4, x^2 - 178, x^3 + 2116]>,
<137, [x + 15, x^2 - 285, x^3 + 5731]>,
<139, [x - 6, x^2 - 176, x^3 - 3008]>,
<149, [x - 2, x^2 - 150, x^3 - 1012]>,
<151, [x - 4, x^2 - 54, x^3 - 342]>,
<157, [x + 3, x^2 - 141, x^3 + 1273]>,
<163, [x - 6, x^2 - 94, x^3 - 1210]>,
<167, [x, x^2 - 262, x^3 + 250]>,
<173, [x, x^2 - 70, x^3 + 406]>,
<179, [x - 3, x^2 - 91, x^3 - 1907]>,
<181, [x + 1, x^2 - 195, x^3 + 1575]>,
<191, [x + 5, x^2 - 213, x^3 + 2729]>,
<193, [x - 12, x^2 - 196, x^3 - 2920]>,
<197, [x + 10, x^2 - 214, x^3 + 4674]>,
<199, [x + 2, x^2 - 86, x^3 + 122]>];

SetPath("QCMod");
load "qc_modular.m";

//z := 1;
//P2<x,y,z> := ProjectiveSpace(Rationals(),2);

/*Q := x^5 - 7/4*x^4*y + 1/8*x^3*y^2 + 7/8*x^2*y^3 - 1/8*x*y^4 - 1/8*y^5 + 13/8*x^3*y*z
    - 17/8*x^2*y^2*z - 1/8*x*y^3*z + 5/8*y^4*z - 1/4*x^3*z^2 + 1/2*x^2*y*z^2 +
    7/8*x*y^2*z^2 - y^3*z^2 + 1/4*x^2*z^3 - 5/8*x*y*z^3 + 5/8*y^2*z^3 +
    1/8*x*z^4 - 1/8*y*z^4;*/
/*Q := x^3*y^2 + 1/2*x^2*y^3 + 1/2*x*y^4 + 1/4*y^5 + 7/4*x^3*y*z + 1/4*x^2*y^2*z +
    5/2*x*y^3*z + 3/2*y^4*z + 1/2*x^3*z^2 - 3/2*x^2*y*z^2 + 3*x*y^2*z^2 +
    5/4*y^3*z^2 - 1/2*x^2*z^3 + 3/2*x*y*z^3 - 1/2*y^2*z^3 + 1/4*x*z^4 -
    1/4*y*z^4;*/
/*Q := x^3*y^2 - 3*x^2*y^3 + 3*x*y^4 - y^5 + 4*x^3*y*z - 8*x^2*y^2*z + 2*x*y^3*z +
    2*y^4*z + 5*x^3*z^2 - x^2*y*z^2 - 5*x*y^2*z^2 + 3*x^2*z^3 - y^2*z^3 + x*z^4;*/
/*Q := x^5 - 1/2*x^4*y + 1/2*x^3*y^2 - 1/2*x^2*y^3 + x^4*z + 3/2*x^3*y*z -
    1/2*x^2*y^2*z + 1/2*x*y^3*z - 1/2*y^4*z - 1/2*x^3*z^2 + 5/2*x^2*y*z^2 +
    1/2*x*y^2*z^2 + y^2*z^3 + 1/2*x*z^4 - 1/2*y*z^4;*/
/*Q := x^5 - 11/2*x^4*y + 97/8*x^3*y^2 - 107/8*x^2*y^3 + 59/8*x*y^4 - 13/8*y^5 -
    2*x^4*z + 6*x^3*y*z - 11/2*x^2*y^2*z + x*y^3*z + 1/2*y^4*z + 9/8*x^3*z^2 -
    9/8*x^2*y*z^2 + 5/8*x*y^2*z^2 - 1/2*y^3*z^2 - 1/8*x^2*z^3 + 1/4*x*y*z^3 +
    1/8*y^2*z^3 - 1/8*x*z^4;*/
/*
X := Scheme(P2, Q);

function TranslationToLineAtInfinity(projective_plane, points)
  if #points notin {2} then error("The number of points must be equal to 2."); end if;
  denom := 1;
  matrix_data := [Eltseq(points[2]), Eltseq(points[1])];
  for x in matrix_data[1] do
    denom := Lcm(Denominator(x), denom);
  end for;
  for x in matrix_data[2] do
    denom := Lcm(Denominator(x), denom);
  end for;
  matrix_data := [[Integers()!(denom * x) : x in elt] : elt in matrix_data];
  points_matrix := Matrix(Rationals(), #points,3, matrix_data);
  // M is a unimodular matrix such that M * P = (1 : 0 : 0) and M * Q = (0 : 1 : 0) lie on the line {z = 0} at infinity
  _, _, M := SmithForm(points_matrix);
  //print M;
  aut := Automorphism(projective_plane, M);
  //print [aut(p) : p in points];
  return aut;
end function;*/

/*
bad_pts := SingularPoints(X);
finite_bad_points := [pt : pt in bad_pts | pt[3] ne 0];
if #finite_bad_points lt 2 then 
    Append(~finite_bad_points, [pt : pt in bad_pts | pt notin finite_bad_points][1]);
end if;
printf "finite bad points: %o\n", finite_bad_points;
aut := TranslationToLineAtInfinity(P2, finite_bad_points);
X2 := aut(X);
print X2;
printf "bad points on X2: %o\n", SingularPoints(X2);
Q := DefiningEquations(X2)[1];
print Q;
exit;
*/
/*Q := x^3*y^2 - 2*x^2*y^3 + x*y^4 + 3*x^3*y*z - 3*x^2*y^2*z - x*y^3*z + y^4*z + 
    5*x^3*z^2 - x^2*y*z^2 - 5*x*y^2*z^2 + 4*x^2*z^3 - 4*x*y*z^3 - 2*y^2*z^3 + 
    x*z^4 - y*z^4;*/ // fourth model
/*Q := x^3*y^2 - x^2*y^3 - x*y^4 + y^5 - 9*x^3*y*z + 17*x^2*y^2*z - 7*x*y^3*z - y^4*z +
    10*x^3*z^2 - 26*x^2*y*z^2 + 19*x*y^2*z^2 - 2*y^3*z^2 + 6*x^2*z^3 - 
    11*x*y*z^3 + 3*y^2*z^3 + x*z^4 - y*z^4; // first model
image_curve_non_monic_eq_xy := Evaluate(Q, [x,y,1]);
print image_curve_non_monic_eq_xy;
exit;*/

/*image_curve_non_monic_eq_xy := x^3*y^2 - 2*x^2*y^3 + x*y^4 + 3*x^3*y - 3*x^2*y^2 - x*y^3 + y^4 + 5*x^3 - x^2*y 
    - 5*x*y^2 + 4*x^2 - 4*x*y - 2*y^2 + x - y; // fourth model*/
/*image_curve_non_monic_eq_xy := x^3*y^2 - x^2*y^3 - x*y^4 + y^5 - 9*x^3*y + 17*x^2*y^2 - 7*x*y^3 
    - y^4 + 10*x^3 - 26*x^2*y + 19*x*y^2 - 2*y^3 + 6*x^2 - 11*x*y + 3*y^2 + x - y; // first model
image_curve_non_monic_eq_xy := image_curve_non_monic_eq_xy / LeadingCoefficient(LeadingCoefficient(image_curve_non_monic_eq_xy));
image_curve_non_monic_eq_xy_lc := LeadingCoefficient(image_curve_non_monic_eq_xy);
image_curve_monic_eq := Numerator(Evaluate(image_curve_non_monic_eq_xy, y / image_curve_non_monic_eq_xy_lc));
print image_curve_monic_eq;
exit;*/

//Q := y^4 + (-2*x^2 - x)*y^3 + (x^4 - 2*x^3 - 8*x^2 - 7*x - 2)*y^2 + (3*x^5 + 5*x^4 - 3*x^3 - 10*x^2 - 6*x - 1)*y + 5*x^6 + 19*x^5 + 28*x^4 + 20*x^3 + 7*x^2 + x; // fourth model
Q := y^5 + (-x - 1)*y^4 + (-x^2 - 7*x - 2)*y^3 + (x^3 + 17*x^2 + 19*x + 3)*y^2 + (-9*x^3 - 26*x^2 - 11*x - 1)*y + 10*x^3 + 6*x^2 + x; // first model

for p_polys in use_polys do
    p, polys := Explode(p_polys);
    /*if p le 7 then
        continue;
    end if;*/
    printf "p = %o\n", p;
    try 
        good_pts, bool, bad_pts, data, fake_rat_pts, bad_disks := QCModAffine(Q, p : use_polys := polys, number_of_correspondences := #polys - 1, N:=20, printlevel:=3);
        printf "good_pts: %o\n", good_pts;
        printf "bool: %o\n", bool;
        printf "bad_pts: %o\n", bad_pts;
        printf "fake_rat_pts: %o\n", fake_rat_pts;
        printf "bad_disks: %o\n\n", bad_disks;
        break;
    catch e 
        print e;
    end try;
end for;