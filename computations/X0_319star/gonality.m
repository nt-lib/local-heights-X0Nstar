X := X0Nstar(319);
//X := X0Nstar(137);
R<x,y,z,w> := Ambient(X);

gon,f := Genus4GonalMap(X);
print f;
f1 := FunctionField(X) ! (DefiningEquations(f)[1]/DefiningEquations(f)[2]);
divs, mult := Support(Divisor(f1));
D := &+[mult[i]*divs[i] : i in [1..#divs] | mult[i] ge 0 ];
H, toX := RiemannRochSpace(CanonicalDivisor(X) -D);
g := toX(H.1)/toX(H.2);
print(g);

P1 := Curve(ProjectiveSpace(Rationals(),1),[]);
f := map< X -> P1 | [x-y+z,y]>;
g := map< X -> P1 | [y,w]>;
assert Degree(f) eq 3;
assert Degree(g) eq 3;

P2<r,s,t> := ProjectiveSpace(Rationals(),2);
pi := map< X -> P2 | [(x+z),w,y]>;
X2 := Image(pi);
F := Evaluate(5*DefiningPolynomial(X2),[y,x,1]); F;