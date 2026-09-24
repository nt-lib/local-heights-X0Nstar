///Magma code accompanying the proof of two results:
// Lemma: If N is squarefree and Phi_2(N) <= 0, then omega(N) <= 8; moreover N<5.44*10^7.
// Theorem: N_0 = 1 231 230=2*3*5*7*11*13*41$. Then N_0 is
// the largest squarefree level with Phi_2(N) <= 0. In particular
// g_0^*(N)>B(N)+1 for every squarefree N>N_0.


// Lemma.
// Theta: differentiation with respect to x.

// a, b, e, C0 are constants.
C<a,b,e,C0> := RationalFunctionField(Rationals(), 4);
F<x> := RationalDifferentialField(C);

// Adjoin sqrt(x), then a formal log(x).
H<Z> := PolynomialRing(F);
K<sqrtx> := ext<F | Z^2-x>;
E<logx> := LogarithmicFieldExtension(K, K!(1/x));

// Work in the final field.
x := E!x;
sqrtx := E!sqrtx;

assert sqrtx^2 eq x;
assert Derivative(x) eq 1;
assert Derivative(sqrtx) eq 1/(2*sqrtx);
assert Derivative(logx) eq 1/x;

theta := a*x - b*sqrtx*(2*logx+C0) - e;

// With the positive square root and the real logarithm:
// Theta(1) = a-b*C0-e < 0
// Since a <= 1/12, e >= 5+13/12, b > 0, C0 > 0.
// Also Theta(x) -> infinity since a > 0.

assert Derivative(theta)
    eq a - b*(2*logx+C0+4)/(2*sqrtx);

assert Derivative(Derivative(theta))
    eq b*(2*logx+C0)/(4*x*sqrtx);

// Thus Theta''(x) > 0 for x >= 1, b > 0, C0 > 0.


// gamma: differentiation with respect to p.

// pb, pi, c0, T, A are constants, as are sqrtP = Sqrt(P) and logP = Log(P).
CP<pb,pi,c0,T,A,sqrtP,logP> := RationalFunctionField(Rationals(), 7);
FP<p> := RationalDifferentialField(CP);

// Adjoin sqrt(p), then a formal log(p).
HP<W> := PolynomialRing(FP);
KP<sqrtp> := ext<FP | W^2-p>;
EP<logp> := LogarithmicFieldExtension(KP, KP!(1/p));

// Work in the final field.
p := EP!p;
sqrtp := EP!sqrtp;
sqrtP := EP!sqrtP;
logP := EP!logP;
P := sqrtP^2;

assert Derivative(sqrtp) eq 1/(2*sqrtp);
assert Derivative(logp) eq 1/p;
assert Derivative(sqrtP) eq 0 and Derivative(logP) eq 0;

// three parts of Phi (at k-th primorial P = P_k), with T = 2^k
Lk := func< P, T | P/(12*T) >;
Mk := func< pb, sqrtP, logP | 2*pb/(4*pi)*sqrtP*(2*logP + c0) >;
Ek := func< k | 5*k + 13/12 >;

// how they change when k -> k+1: P -> p*P, T -> 2*T, pb -> (1/2 + 1/sqrt(p))*pb
assert Lk(p*P, 2*T)/Lk(P, T) eq p/2;

gamma := func< sqrtp, logp, A | (1/sqrtp + 2/sqrtp^2)*(1 + 2*logp/A) >;
assert Mk((1/2 + 1/sqrtp)*pb, sqrtp*sqrtP, logP + logp)/(Mk(pb, sqrtP, logP)*p/2)
    eq gamma(sqrtp, logp, 2*logP + c0);

// derivative of gamma w.r.t. to p is as claimed, so gamma is indeed decreasing in p
assert Derivative(gamma(sqrtp, logp, A))/gamma(sqrtp, logp, A)
    eq (1/2)/(p + 2*sqrtp) + 2/(p*(A + 2*logp)) - 1/p;

// direct check from the line above it is < 0
// for p >= 1 (so logp >= 0) and A > 4, since then 2/(A+2*logp) < 1/2


// the next is just a confirmation that three parts of Phi were correctly defined
RR := RealField(30);
C0 := 10 - 4*Log(RR!3);
pb := func< k | &*[RR| 1/2 + 1/Sqrt(RR!NthPrime(i)) : i in [1..k]] >;
Pk := func< k | &*[Integers()| NthPrime(i) : i in [1..k]] >;
Ak := func< k | 2*Log(RR!Pk(k)) + C0 >;
Phi := func< x, k | (RR!x)/(12*2^k) - 2*pb(k)/(4*Pi(RR))*Sqrt(RR!x)*(2*Log(RR!x) + C0) - 5*k - (RR!13)/12 >;
Lk := func< P, k | (RR!P)/(12*2^k) >;
Mk := func< pb, P | 2*pb/(4*Pi(RR))*Sqrt(RR!P)*(2*Log(RR!P) + C0) >;
Ek := func< k | 5*k + (RR!13)/12 >;
Abs(Phi(Pk(9), 9) - (Lk(Pk(9), 9) - Mk(pb(9), Pk(9)) - Ek(9))) lt (RR!1)/10^20;
// -> true

gamma := func< p, A | (1/Sqrt(RR!p) + 2/(RR!p))*(1 + 2*Log(RR!p)/A) >;
gamma(NthPrime(10), Ak(9));
// -> 0.293593178000983694262562428554

// base case
Phi(Pk(9), 9);
// -> 9186.32678401491177834604328637


// Theorem.
// we now move to computing the largest squarefree level N with Phi_2(N) negative
// Phi_2(N) for squarefree N, from its factorisation Fac = Factorisation(N)
PhiF := function(N, Fac)
    psi := &*[Integers()| t[1]+1 : t in Fac];
    pi_k := &*[RR| 1/2 + 1/Sqrt(RR!t[1]) : t in Fac];
    return psi/(12*2^#Fac) - 2*pi_k/(4*Pi(RR))*Sqrt(RR!N)*(2*Log(RR!N) + C0) - 5*#Fac - (RR!13)/12;
end function;

PhiN := func< N | PhiF(N, Factorisation(N)) >;

// Z_{2,k}, the unique zero of Phi_{2,k}, by bisection/binary search
// Magma doesn't seem to have a function supporting this!
Zk := function(k)
    lo := RR!1; hi := RR!10^13;
    for i in [1..300] do
        m := (lo + hi)/2;
        if Phi(m, k) gt 0 then hi := m; else lo := m; end if;
    end for;
    return hi;
end function;

// Phi_2(N) <= 0 forces omega(N) <= 8 and N <= Z_{2,omega(N)}, and Z_{2,k} increases in k,
// so nothing above Z_{2,8} can qualify 
Bound := Ceiling(Zk(8));
Bound;
// -> 54373185

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
// -> 1231230 9605

PhiN(1231230);
// -> -17.8149833030895690414382846436
