
QQ := Rationals();
R<x, y> := PolynomialRing(QQ, 2);
Q := y^4 + (2*x^6 - 53/6*x^5 + 179/18*x^4 - 17/9*x^3 - 53/18*x^2 + 16/9*x -
5/18)*y^3 + (1/2*x^12 - 15/4*x^11 + 547/324*x^10
    + 3826/81*x^9 - 28171/162*x^8 + 50029/162*x^7 - 12079/36*x^6 + 78119/324*x^5
- 19073/162*x^4 + 3151/81*x^3 -
    1357/162*x^2 + 173/162*x - 5/81)*y^2 + (13/18*x^17 - 2215/162*x^16 +
629557/5832*x^15 - 461353/972*x^14 +
    960544/729*x^13 - 177451/72*x^12 + 18654127/5832*x^11 - 695638/243*x^10 +
1190911/729*x^9 - 2121427/5832*x^8 -
    220879/729*x^7 + 92491/243*x^6 - 54091/243*x^5 + 61304/729*x^4 -
31405/1458*x^3 + 886/243*x^2 - 10/27*x +
    25/1458)*y + 5/18*x^22 - 809/108*x^21 + 173117/1944*x^20 -
65354405/104976*x^19 + 25284455/8748*x^18 -
    1003358513/104976*x^17 + 15274865/648*x^16 - 1566356959/34992*x^15 +
439320299/6561*x^14 - 8418112045/104976*x^13 +
    170060350/2187*x^12 - 358620545/5832*x^11 + 1044449483/26244*x^10 -
368832007/17496*x^9 + 119482381/13122*x^8 -
    83483947/26244*x^7 + 1938787/2187*x^6 - 1686341/8748*x^5 + 207152/6561*x^4 -
16055/4374*x^3 + 1775/6561*x^2 -
    125/13122*x;

p := 29;
//reduce the curve mod p
Fp := GF(p);
R_p<x_p, y_p> := PolynomialRing(Fp, 2);
phi := hom<R -> R_p | [x_p, y_p]>;
Q_p := phi(Q);
print "The reduction of the curve modulo p is Q_p = 0, where Q_p is:";
print Q_p;

dQdx := Derivative(Q_p, x_p);
dQdy := Derivative(Q_p, y_p);

rational_points := PointSearch(Curve(AffineSpace(QQ, 2), Q), 1000);
rational_points;
for P in rational_points do
    px := P[1];
    py := P[2];
    // Reduce the point's coordinates modulo p
    px_p := Fp ! px;
    py_p := Fp ! py;
    P_p := [px_p, py_p];
    printf "Rational point P = (%o, %o) reduces to P_p = (%o, %o)\n", px, py,
px_p, py_p;
    // Sanity check: does the reduced point lie on the reduced curve?
    if Evaluate(Q_p, P_p) ne 0 then
        printf "  --> WARNING: The reduced point does NOT lie on the reduced
curve. Check your point coordinates.\n";
        continue;
    end if;
    // Evaluate partial derivatives at the reduced point
    val_dx := Evaluate(dQdx, P_p);
    val_dy := Evaluate(dQdy, P_p);
    printf "  Partial derivatives at P_p: dQ/dx = %o, dQ/dy = %o\n", val_dx, val_dy;
    // Check the smoothness condition
    if val_dx eq 0 and val_dy eq 0 then
        printf "  --> RESULT: The reduced point is SINGULAR.\n\n";
    else
        printf "  --> RESULT: The reduced point is SMOOTH.\n\n";
    end if;
end for;
