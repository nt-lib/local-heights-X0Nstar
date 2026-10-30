// NOTE: the output of this file is not deterministic.  It chooses the linear forms at random, so
// every run finds different models (or none).  The plane model used in qc_p37_mws.m was found by
// such a random search; its linear forms L1, L2, L3, L4 are hardcoded there, so the proof does not
// depend on this file.
//
// This file is an example of a random search for plane models of X_0(187)^*: it maps the canonical
// model of X_0(187)^* to plane models given by random linear forms with coefficients in [-2..2] and
// reports those that are good for quadratic Chabauty at p = 7 (no F_7-points of the canonical model
// map to bad or infinite points of the plane model).  It stops after finding five such models or
// after 100000 trials.  See random_model_search in src/model_search.m.

SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
load "src/model_search.m";

X187star := X0Nstar(187);
pts := PointSearch(X187star, 30);
// name the coordinates, so that the linear forms found are printed in X, Y, Z
P<X, Y, Z> := AmbientSpace(X187star);

p := 7;
models := random_model_search(X187star, pts, p, x, y : coefficient_bound := 2, max_trials := 100000,
  max_successes := 5);

if #models eq 0 then
  printf "No model that is good at p = %o found.\n", p;
else
  printf "Found %o models that are good at p = %o.\n", #models, p;
end if;
