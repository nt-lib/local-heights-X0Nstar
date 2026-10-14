AttachSpec("local_heights.spec");

Ns := [161, 187, 247, 319, 290];
for N in Ns do
    printf "N = %o:\n", N;
    NTLStarQuotientDualGraph(N);
    printf "\n";
end for;
