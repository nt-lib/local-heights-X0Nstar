/*
Random search for plane models that are good for quadratic Chabauty at a given prime.
Not part of the local_heights.spec package; load it after ModularCurvesX0plusG4-6/model_equation_finder.m.
*/

// Randomly search for plane models of domain_curve (a canonical model in P^(g-1) with known rational
// points known_pts) given by linear forms L1, L2, L3, L4 in the coordinates of its ambient space, with
// integer coefficients in [-coefficient_bound..coefficient_bound]; the model is the image of the map
// (L1/L2, L3/L4) as in find_and_test_model. A model is good if p is a good prime for it in the sense of
// find_and_test_model: no F_p-points of domain_curve map to bad or infinite points of the model.
// x and y are the variables for the model equation, as in find_and_test_model.
// Tries at most max_trials random choices and stops after finding max_successes good models.
// The linear forms are chosen at random, so the result differs between runs.
// Returns the sequence of the good models found, each as the sequence [L1, L2, L3, L4].
function random_model_search(domain_curve, known_pts, p, x, y : coefficient_bound := 2, max_trials := 100000,
    max_successes := 5, printlevel := 1)
  R := CoordinateRing(AmbientSpace(domain_curve));
  coeffs := [-coefficient_bound..coefficient_bound];
  lins := [&+[c[i]*R.i : i in [1..Rank(R)]] : c in CartesianPower(coeffs, Rank(R))];
  lins := [L : L in lins | L ne 0];

  found := [];
  for trial in [1..max_trials] do
    Ls := [Random(lins) : j in [1..4]];
    // x or y would be constant
    if Ls[1] eq Ls[2] or Ls[3] eq Ls[4] then
      continue;
    end if;
    try
      _, _, _, _, good_primes := find_and_test_model(Ls[1], Ls[2], Ls[3], Ls[4], domain_curve, x, y, known_pts :
        min_prime := p, max_prime := p, printlevel := printlevel);
    catch e
      // e.g. the map to P^2 is not birational onto its image
      continue;
    end try;
    if p in good_primes then
      if printlevel ge 1 then
        printf "Trial %o: found a model that is good at p = %o: L1 = %o, L2 = %o, L3 = %o, L4 = %o\n",
          trial, p, Ls[1], Ls[2], Ls[3], Ls[4];
      end if;
      Append(~found, Ls);
      if #found ge max_successes then
        break;
      end if;
    end if;
  end for;
  return found;
end function;
