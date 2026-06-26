load "src/StarQuotientMeasures.m";

Ns := [161, 187, 247, 319];
for N in Ns do
    printf "N = %o:\n", N;
    StarQuotientDualGraph(N, [2]);
    printf "\n";
end for;