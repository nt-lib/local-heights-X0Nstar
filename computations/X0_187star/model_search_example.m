// NOTE: the output of this file is not deterministic.  It chooses the linear forms at random, so
// every run finds different models (or none).  This file is only meant to show how to find new
// models, not to reproduce ours.  The plane model used in qc_p37_mws.m was found once by a similar
// random search (not targeting a specific prime), and the linear forms L1, L2, L3, L4 found are
// hardcoded in qc_p37_mws.m, so that our computations are reproducible.
//
// This file is an example of a random search for plane models of X_0(187)^*: it maps the canonical
// model of X_0(187)^* to plane models given by random linear forms and reports those that have a
// good prime for quadratic Chabauty (a prime p such that no F_p-points of the canonical model map
// to bad or infinite points of the plane model).  See random_model_search in src/model_search.m
// for the default search parameters.

SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
load "src/model_search.m";

X187star := X0Nstar(187);
pts := PointSearch(X187star, 30);
// name the coordinates, so that the linear forms found are printed in X, Y, Z
P<X, Y, Z> := AmbientSpace(X187star);

models := random_model_search(X187star, pts, x, y);
printf "Found %o models with a good prime.\n", #models;
