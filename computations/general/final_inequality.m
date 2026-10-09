///Magma code accompanying the proof of two results:
// Lemma: If N is squarefree and Phi(N) <= 0, then omega(N) <= 9; moreover N<2.87*10^8.
// Theorem: N_0 = 2 589 510=2*3*5*7*11*19*59. Then N_0 is
// the largest squarefree level with Phi(N) <= 0. In particular
// g_0^*(N)>B(N)+1 for every squarefree N>N_0.
// N_0 must equal N_max in src/genera_and_bounds.m, up to which existence_of_relations.m
// checks the remaining levels directly; this is asserted below.

load "src/genera_and_bounds.m";   // for N_max

// Lemma.
RR := RealField(30);
C0 := 10 - 4*Log(RR!3);
assert C0 gt 4;
pb := func< k | &*[RR| 1/2 + 1/Sqrt(RR!NthPrime(i)) : i in [1..k]] >;
Pk := func< k | &*[Integers()| NthPrime(i) : i in [1..k]] >;
Ak := func< k | 2*Log(RR!Pk(k)) + C0 >;
Phi := func< x, k | (RR!x)/(12*2^k) - 3*pb(k)/(4*Pi(RR))*Sqrt(RR!x)*(2*Log(RR!x) + C0) - 5*k - (RR!13)/12 >;

gamma := func< p, A | (1/Sqrt(RR!p) + 2/(RR!p))*(1 + 2*Log(RR!p)/A) >;
gamma(NthPrime(10), Ak(9));
// -> 0.293593178000983694262562428554
assert gamma(NthPrime(10), Ak(9)) lt 1;

// base case
Phi(Pk(10), 10);
// -> 353541.447464298701224548679871
assert Phi(Pk(10), 10) gt 0;

// Z_k, the unique zero of Phi_k, by bisection/binary search
// Magma doesn't seem to have a function supporting this!
Zk := function(k)
    lo := RR!1; hi := RR!10^13;
    assert Phi(lo, k) lt 0 and Phi(hi, k) gt 0;
    for i in [1..300] do
        m := (lo + hi)/2;
        if Phi(m, k) gt 0 then hi := m; else lo := m; end if;
    end for;
    return hi;
end function;

// Z_1 < ... < Z_9 < 2.87*10^8
Zs := [Zk(k) : k in [1..9]];
assert &and[Zs[k] lt Zs[k+1] : k in [1..8]] and Zs[9] lt 287*10^6;


// Theorem.
// we now move to computing the largest squarefree level N with Phi(N) negative
// Phi(N) for squarefree N, from its factorisation Fac = Factorisation(N).
// The paper only defines Phi(N) for N > 6.
PhiF := function(N, Fac)
    error if N le 6, Sprintf("Phi(N) is only defined for N > 6, not for N = %o.", N);
    psi := &*[Integers()| t[1]+1 : t in Fac];
    pi_k := &*[RR| 1/2 + 1/Sqrt(RR!t[1]) : t in Fac];
    return psi/(12*2^#Fac) - 3*pi_k/(4*Pi(RR))*Sqrt(RR!N)*(2*Log(RR!N) + C0) - 5*#Fac - (RR!13)/12;
end function;

PhiN := func< N | PhiF(N, Factorisation(N)) >;

// Phi(N) <= 0 forces omega(N) <= 9 and N <= Z_{omega(N)}, and Z_k increases in k,
// so nothing above Z_9 can qualify
Bound := Ceiling(Zs[9]);
Bound;
// -> 286055291

// test every N up to the bound
largest := 0; cnt := 0;
for N in [7..Bound] do
    if N mod 10000000 eq 0 then printf "  %o  %o s\n", N, Cputime(); end if;
    if N mod 4 eq 0 or N mod 9 eq 0 or N mod 25 eq 0 then continue; end if;
    Fac := Factorisation(N);
    if forall{t : t in Fac | t[2] eq 1} and PhiF(N, Fac) le 0 then
        cnt +:= 1;
        largest := N;
    end if;
end for;
largest, cnt;
// -> 2589510 23189
assert largest eq N_max;

PhiN(N_max);
// -> -15.0507895145218108862184033655
assert PhiN(N_max) le 0;


// Facts used in the proof of the lemma, easily checked by hand:
// d/dx ( sqrt(x)*(2*log(x) + C0) )           = (2*log(x) + C0 + 4)/(2*sqrt(x))
// d/dx ( (2*log(x) + C0 + 4)/(2*sqrt(x)) )   = -(2*log(x) + C0)/(4*x^(3/2)) < 0 for x >= 1,
//     so Theta' is increasing; Theta(1) = a - b*C0 - e < 0 and Theta(x) -> oo as x -> oo.
// With L_k = P_k/(12*2^k), M_k = 3*pb_k/(4*pi)*sqrt(P_k)*(2*log(P_k) + C0), A_k = 2*log(P_k) + C0:
//     L_{k+1} = (p/2)*L_k,  M_{k+1} = gamma(p, A_k)*(p/2)*M_k,  A_{k+1} = A_k + 2*log(p).
// d/dp log(gamma(p, A)) = (1/2)/(p + 2*sqrt(p)) + 2/(p*(A + 2*log(p))) - 1/p
//     <= (1/p)*(2/A - 1/2) < 0 for p >= 1, A > 4, so gamma is decreasing in p (and in A).


// For the reader who still wants this confirmed in Magma, uncomment the following:
/*
// Theta: differentiation with respect to x, with a, b, e, C0 constants;
// sqrt(x) and log(x) are adjoined as an algebraic and a logarithmic extension.
C<a,b,e,C0> := RationalFunctionField(Rationals(), 4);
F<x> := RationalDifferentialField(C);
H<Z> := PolynomialRing(F);
K<sqrtx> := ext<F | Z^2-x>;
E<logx> := LogarithmicFieldExtension(K, K!(1/x));
x := E!x;
sqrtx := E!sqrtx;
assert sqrtx^2 eq x;
assert Derivative(x) eq 1;
assert Derivative(sqrtx) eq 1/(2*sqrtx);
assert Derivative(logx) eq 1/x;

theta := a*x - b*sqrtx*(2*logx+C0) - e;
assert Derivative(sqrtx*(2*logx+C0)) eq (2*logx+C0+4)/(2*sqrtx);
assert Derivative((2*logx+C0+4)/(2*sqrtx)) eq -(2*logx+C0)/(4*x*sqrtx);
assert Derivative(theta) eq a - b*(2*logx+C0+4)/(2*sqrtx);
assert Derivative(Derivative(theta)) eq b*(2*logx+C0)/(4*x*sqrtx);

// gamma: differentiation with respect to p. As P = P_k is constant in p, sqrt(P) and log(P)
// are separate constants sqrtP, logP (with P = sqrtP^2); T = 2^k.
CP<pb,pi,c0,T,A,sqrtP,logP> := RationalFunctionField(Rationals(), 7);
FP<p> := RationalDifferentialField(CP);
HP<W> := PolynomialRing(FP);
KP<sqrtp> := ext<FP | W^2-p>;
EP<logp> := LogarithmicFieldExtension(KP, KP!(1/p));
p := EP!p;
sqrtp := EP!sqrtp;
sqrtP := EP!sqrtP;
logP := EP!logP;
P := sqrtP^2;
assert Derivative(sqrtp) eq 1/(2*sqrtp);
assert Derivative(logp) eq 1/p;
assert Derivative(sqrtP) eq 0 and Derivative(logP) eq 0;

// L_k, M_k and how they change for k -> k+1:
// P -> p*P, T -> 2*T, pb -> (1/2 + 1/sqrt(p))*pb, log(P) -> log(P) + log(p)
Lk := func< P, T | P/(12*T) >;
Mk := func< pb, sqrtP, logP | 3*pb/(4*pi)*sqrtP*(2*logP + c0) >;
assert Lk(p*P, 2*T)/Lk(P, T) eq p/2;
gamma := func< sqrtp, logp, A | (1/sqrtp + 2/sqrtp^2)*(1 + 2*logp/A) >;
assert Mk((1/2 + 1/sqrtp)*pb, sqrtp*sqrtP, logP + logp)/(Mk(pb, sqrtP, logP)*p/2)
    eq gamma(sqrtp, logp, 2*logP + c0);
assert Derivative(gamma(sqrtp, logp, A))/gamma(sqrtp, logp, A)
    eq (1/2)/(p + 2*sqrtp) + 2/(p*(A + 2*logp)) - 1/p;
print "All derivative identities confirmed.";
*/
