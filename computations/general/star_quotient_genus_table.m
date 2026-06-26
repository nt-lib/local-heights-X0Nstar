load "src/StarQuotientMeasures.m";


data := [];
for N in [2..999] do
    if not IsSquarefree(N) then continue; end if;
    g := StarQuotientGenus(N);
    c := NrOfStarQuotientComponentConditions(N);
    Append(~data, <g, N, c>);
end for;

Sort(~data);

printf "N\tgenus\tnr_conditions\n";
for t in data do
    printf "%o\t%o\t%o\n", t[2], t[1], t[3];
end for;

for t in data do
    g, N, c := Explode(t);
    if g - c ge 2 then
        continue;
    end if;
    r_list := [p : p in PrimesUpTo(100) | not IsDivisibleBy(N, p)];
    polys, is_hecke_generator := StarQuotientGoodCorrespondences(N, r_list);
    print polys, is_hecke_generator; // TODO: empty polys for N = 249
    if not exists{i : i in [1..#polys] | is_hecke_generator[i]} then
        printf "For N = %o, there is no poly in T_r.\n", N;
    else
        printf "N = %o violates the inequality, but there exists a non-trivial poly.\n", N;
    end if;
end for;