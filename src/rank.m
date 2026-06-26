// This file contains functions for computing the ranks (over Q) of J_0(N) and its factors

// The code contains the following functions (further descriptions before each function)

// simple_rank: Attempts to compute the rank of a simple factor of J_0(N)
// rank_quo: Attempts to compute the rank of an Atkin-Lehner quotient of the Jacobian of J_0(N)
// equal_rank: Attempts to check whether the rank of J_0(N) and a given quotient are equal 

// There are examples at the end of the file

///////////////////
/// simple_rank ///
///////////////////

// Input : A
// A is a simple modular abelian variety that appears as a factor of J_0(N)

// Ouput : Rank(A), tf
// Rank(A) and whether this is provably true.
// If the second output is false then Rank(A) will be -1000

simple_rank := function(A);
    tf := true;
    if Dimension(A) eq 0 or IsZeroAt(LSeries(A),1) eq false then
       rk := 0;
    else M := ModularSymbols(Newform(A));
         _,ord_vanishing := LSeriesLeadingCoefficient(M,1,100);
         if ord_vanishing eq 1 then
            rk := Degree(HeckeEigenvalueField(M));
         else if Dimension(A) eq 1 then   // we have an elliptic curve and we try to compute the rank directly
                 E := EllipticCurve(A);
                 rk,ell_tf := Rank(E);
                 assert ell_tf eq true;
              else rk := -1000;
                   tf := false;
              end if;
          end if;
    end if;
    return rk, tf;
end function;

////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////
/// rank_quo ///
////////////////

// Input: N, sequence of AL indices which generates a group W
// Ouput: Rank(J_0(N)) / W, tf
// Rank and whether it is provably true
// If the second output is false then the rank will be -1000

// if you want to compute the rank of J_0(N) (without quotienting) then take empty sequence

rank_quo := function(N,seq_al);
    J := JZero(N);
    if #seq_al gt 0 then
        J_quo := Image(&*[1+AtkinLehnerOperator(J,i) : i in seq_al]);
    else J_quo := J;
    end if;
    dec := Decomposition(J_quo);
    rk := 0;
    tf := true;

    for A in dec do
        rkA, Atf := simple_rank(A);
        if Atf eq false then  // cannot compute rank with simple_rank
           rk := -1000;
           tf := false;
           break;
        else rk := rk + rkA;
        end if;
     end for;

     return rk, tf;
end function;