//this code is performing a random/arbitrary search
//it maps a canonical model of X0(187)^* into
//a bunch of plane affine models
// and checks which of them are good for QC prime 7
// (good in the sense that there are no F_7-points
//   on the canonical model mapping into bad
//   or infinite disks)
// It stops at most maxTrials (=100 000) attempts
// or after finding five such models.


// N. B. One also needs to do a two-line modification
// of model_equation_finder.m to check only prime 7
// (otherwise it checks all primes up to 50)
// Due to randomness, the output changes with re-runs.

SetPath("QCMod");
load "ModularCurvesX0plusG4-6/model_equation_finder.m";


X187star := X0Nstar(187);
pts := PointSearch(X187star, 30);

// Ambient projective space
P<X, Y, Z> := ProjectiveSpace(Rationals(), 2);

// Function to build all linear forms with coeffs in [-2..N]
LinearForms := function(N)
    forms := [];
    for a in [-2..N] do
    for b in [-2..N] do
    for c in [-2..N] do
        Append(~forms, a*X + b*Y + c*Z);
    end for; end for; end for;
    return forms;
end function;

// Generate search space
N := 2;
lins := LinearForms(N);

success := 0;
count := 0;
maxTrials := 100000;  // number of random trials


for trial in [1..maxTrials] do
    // pick 4 random linear forms
    Ls := [ Random(lins) : j in [1..4] ];
    L1 := Ls[1]; L2 := Ls[2]; L3 := Ls[3]; L4 := Ls[4];

    // apply the same "continue" rules as before
    if L1 eq L2 or L3 eq L4 then
        continue;
    end if;
                count +:= 1;

                if IsPrime(count) and count mod 12 eq 1 then
    printf "Trying combination %o: L1=%o, L2=%o, L3=%o, L4=%o\n", count, L1, L2, L3, L4;
end if;

                try
                    S, Q, model_map, pts187, good_primes, inf_pts, bad_pts :=
                        find_and_test_model(L1, L2, L3, L4, X187star, x, y, pts : max_prime:=7, printlevel:=1);
                catch e
                    continue;
                end try;

                if #good_primes gt 0 then
			success := success +1;
                    print "SUCCESS! Good primes found.";
                    print "L1 =", L1;
                    print "L2 =", L2;
                    print "L3 =", L3;
                    print "L4 =", L4;
                end if;
  	if success gt 5 then break; end if;
end for;

if success eq 0 then
    print "No suitable model found *randomly* up to coefficient bound", N;
end if;
