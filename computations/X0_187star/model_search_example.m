/* 
This file is an example of a random search for plane models of X_0(187)^*: it maps the canonical
model of X_0(187)^* to plane models given by random linear forms and reports those that have a
good prime for quadratic Chabauty (a prime p such that no F_p-points of the canonical model map
to bad or infinite points of the plane model).  See random_model_search in src/model_search.m
for the default search parameters.

NOTE: the output of this file is not deterministic. A random set of linear forms L1, L2, L3 and L4
that produced a particularly nice model mod 37 have been hardcoded in qc_p37_mws.m in order to make
the rational point computations in this repository reproducible. Searches similar to those in this 
file have also been carried out for X_0(247)^*, X_0(319)^* and X_0(430)^*. 
*/

SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";
load "src/model_search.m";

X187star := X0Nstar(187);
pts := PointSearch(X187star, 30);
// name the coordinates, so that the linear forms found are printed in X, Y, Z
P<X, Y, Z> := AmbientSpace(X187star);

models := random_model_search(X187star, pts, x, y);
printf "Found %o models with a good prime.\n", #models;
