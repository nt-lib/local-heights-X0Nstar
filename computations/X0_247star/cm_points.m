// Rational CM points on X_0(247)^*.
// RationalCMDiscs from ShimuraCurveALQuotients uses Shimura reciprocity to find the
// discriminants of the orders whose CM points on X_0(N)^* are rational, and the number of
// rational points for each. These CM points are distinct rational points other than the
// cusp. qc_p31.m proves that the known rational points are all rational points, so if
// the counts below agree, every rational point other than the cusp is a CM point, with one
// of the discriminants below.
AttachSpec("ShimuraCurveALQuotients/ShimuraQuotients.spec");

N := 247;
Xstar := CreateShimuraQuot(1, N, Set(Divisors(N)));
discs := RationalCMDiscs(Xstar);
printf "Discriminants of the rational CM points on X_0(%o)^*, with the number of points:\n%o\n",
    N, [<d, discs[d]> : d in Reverse(Sort(Setseq(Keys(discs))))];

// For squarefree N the curve X_0(N)^* has exactly one cusp.
pts := PointSearch(X0Nstar(N), 30);
printf "Known rational points on X_0(%o)^* (PointSearch, height bound 30):\n%o\n", N, pts;
assert #pts eq 1 + &+[discs[d] : d in Keys(discs)];
printf "The %o known rational points are the cusp and %o CM points.\n", #pts, #pts - 1;
