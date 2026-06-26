SetPath("QCMod");
load "ModularCurvesX0plusG4-6/mef_new_modified.m";

max_prime := 100;

X319star := X0NQuotient(319, [11, 29]);
pts_known := PointSearch(X319star, 100);

printf "Found %o rational points on the curve:\n%o\n\n", #pts_known, pts_known;

// Loop through each rational point, using it as a projection center.
for i in [1..#pts_known] do
    pt := pts_known[i];
    print "========================================================================";
    printf "===== Starting Search %o/%o: Projecting from point %o =====\n", i, #pts_known, pt;
    print "========================================================================";

    // Define the initial map from the curve to a plane model.
    image_C_scheme, projection_map := ProjectionFromNonsingularPoint(X319star, pt);

    // The 'projection_map' returned is a composite map. We need to extract the
    // final defining polynomials to construct a valid map object for our function.
    P2 := AmbientSpace(image_C_scheme);
    C_clean := Curve(P2, DefiningEquation(image_C_scheme));

    map_parts := Components(projection_map);
    assert #map_parts eq 2;
    automorphism_polys := DefiningPolynomials(map_parts[1]);
    final_defining_polys := [ automorphism_polys[j] : j in [2..4] ];

    initial_model_map := map< X319star -> C_clean | final_defining_polys >;



    print "Starting the search for a good plane model...";
    // Assign the multiple return values directly to a sequence of variables.
    image_curve, image_curve_eq, model_map, image_ratpts, good_primes,
    infinite_primes_data, bad_primes_data, bad_disks_eq, _, _
        := find_and_test_model_2(initial_model_map, pts_known : max_prime := max_prime, printlevel := 1);

    // On failure, the function returns empty lists. Check for the success type.
    if Type(image_curve) eq CrvPln then

        print "\n\n-----------------------------------------------";
        printf "          RESULTS FOR PROJECTION FROM %o\n", pt;
        print "-----------------------------------------------\n";

        printf "A suitable model was found.\n\n";
        printf "The new affine equation for the curve is:\n%o = 0\n\n", image_curve_eq;
        printf "The map from the original curve to the new model in P^2 is:\n%o\n\n", model_map;
        printf "The images of the known rational points on the new model are:\n%o\n\n", image_ratpts;
        printf "The polynomial whose roots are the 'bad' x-coordinates is:\n%o\n\n", bad_disks_eq;
        printf "For this model, the following primes up to %o are 'good' (no bad points or points at infinity mod p):\n%o\n", max_prime, good_primes;

        if #infinite_primes_data gt 0 then
            print "\n--- Primes where points map to infinity on the image curve ---";
            for item in infinite_primes_data do
                printf "p = %o: %o points\n", item`p, #item`infinite_points;
            end for;
        end if;

        if #bad_primes_data gt 0 then
            print "\n--- Primes where the reduction has 'bad' points ---";
            for item in bad_primes_data do
                 printf "p = %o: %o points\n", item`p, #item`bad_points;
            end for;
        end if;
    else
        print "\n\n-----------------------------------------------";
        printf "      NO SUITABLE MODEL FOUND FOR PROJECTION FROM %o\n", pt;
        print "-----------------------------------------------";
        if Type(image_curve) ne CrvPln then
            print "The search function exited early";
        else 
             printf "A model was generated, but it has no 'good' primes up to the limit of %o.", max_prime;
        end if;
    end if;
    print "========================================================================";
end for;

