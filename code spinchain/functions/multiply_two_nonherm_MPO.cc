#include "itensor/all.h"
#include "header.h"

using namespace itensor;

// TODO Leon fragen: stimmt das mit dem orthocenter so?
// TODO Leon fragen: hier besser factor verwenden?

MPO liouvillian::multiply_two_nonherm_MPO(const MPO& MPO_A, const MPO& MPO_B)
// Beachte bei der Eingabe: Die beiden MPOs müssen mit den selben Site-Indizes versehen sein. Bei beiden muss der 
// Zeilenindex ungeprimed und der Spaltenindex geprimed sein. 
{
    // set up empty MPS for result
    int N = MPO_A.N();
    int K = MPO_B.N(); 
    if (N!=K) 
    { 
        std::cout << "In function multiply_two_nonherm_MPO: size of MPOs does not fit. Abort" << std::endl; 
    }
    
    MPO A = MPO_A;
    MPO B = MPO_B;
    MPO res = MPO(N);
    
    A.position(1);
    B.position(1);
    
    B.primeall();
    
    Index si_A;
    Index li_res;
    ITensor temp, fork, U, S, V;
    auto IndexTypeSiteNoPrime = [](Index const& i) { return (i.type() == Site && i.primeLevel() == 0); };
    
    // loop through sites. First and last site are special cases
    for (int i = 1; i<=N; i++)
    {
        if (i == 1)
        { 
            si_A = findindex(A.A(i), IndexTypeSiteNoPrime);
            temp = A.A(i) * B.A(i);
            U = ITensor(si_A, prime(si_A,2));
            svd(temp, U, S, V);
            U = U * S;
            res.setA(i,U);
        }
        else if (i == N)
        {
            si_A = findindex(A.A(i), IndexTypeSiteNoPrime);
            li_res = commonIndex(U, V, Link);
            temp = V * A.A(i) * B.A(i);
            res.setA(i,temp);
        }
        else
        {
            si_A = findindex(A.A(i), IndexTypeSiteNoPrime);
            li_res = commonIndex(U, V, Link);
            temp = V * A.A(i) * B.A(i);
            U = ITensor(si_A, prime(si_A,2), li_res);
            svd(temp, U, S, V);
            U = U * S;
            res.setA(i,U);
        }
        
    } // EOF site loop
    
    res.mapprime(2,1,Site);
    res.orthogonalize();
    return res;  
}











