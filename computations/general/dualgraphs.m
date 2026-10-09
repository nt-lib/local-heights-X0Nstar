load "src/StarQuotientMeasures.m";

Ns := [161, 187, 247, 319, 290];
for N in Ns do
    printf "N = %o:\n", N;
    StarQuotientDualGraph(N);
    printf "\n";
end for;
