// The genus of X_0(5001)^*, and the time it takes to compute a basis of its clean
// correspondences as polynomials in T_q.
AttachSpec("local_heights.spec");

N := 5001;
g := NTLGenusX0Nstar(N);
assert GenusX0NQuotient(N, PrimeDivisors(N)) eq g;
printf "Genus of X_0(%o)^*: %o\n", N, g;
assert g eq 117;

// The smallest prime q coprime to N such that T_q generates the Hecke algebra.
S, Sstar := NTLStarQuotientCuspForms(N);
assert exists(q){q : q in PrimesUpTo(100) | not IsDivisibleBy(N, q)
    and NTLIsHeckeGenerator(NTLHeckeMatrixOnFixedSubspace(S, Sstar, q))};
printf "T_%o generates the Hecke algebra.\n", q;

t := Cputime();
polys := NTLStarQuotientCleanCorrespondences(N, [q] : trace_zero := true)[1];
printf "Basis of the clean correspondences as polynomials in T_%o: %o polynomials of degrees %o, computed in %o seconds.\n",
    q, #polys, [Degree(f) : f in polys], Cputime(t);
