#include "itensor/all.h"
#include "header.h"

namespace lv = liouvillian;

double liouvillian::calc_rel_current_new(MPS& rho, const int& m, const std::vector<Index>& s)
{
   // get size of rho
    int N = rho.N();
    
    // initialize current and tensors
    ITensor temp;
    double rel_curr = 0.0;
    
    // initialize
    MPO CurrMPO1 = MPO(N);
    MPO CurrMPO2 = MPO(N);
    ITensor MPOTensor1, MPOTensor2;
    
    // set up operators 
    for (int i=0; i<N; i++)
    {
        if (i == m) 
        { 
            MPOTensor1 = ITensor(s[i],prime(s[i])); 
            MPOTensor1.set(s[i](2),prime(s[i](2)),1.0);
            MPOTensor2 = ITensor(s[i],prime(s[i])); 
            MPOTensor2.set(s[i](3),prime(s[i](3)),1.0);
        }
        else if (i == (m+1)) 
        { 
            MPOTensor1 = ITensor(s[i],prime(s[i])); 
            MPOTensor1.set(s[i](3),prime(s[i](3)),1.0);
            MPOTensor2 = ITensor(s[i],prime(s[i])); 
            MPOTensor2.set(s[i](2),prime(s[i](2)),1.0);
        }
        else 
        { 
            MPOTensor1 = ITensor(s[i],prime(s[i]));
            MPOTensor1.set(s[i](1),prime(s[i](1)),1.0);
            MPOTensor2 = ITensor(s[i],prime(s[i])); 
            MPOTensor2.set(s[i](1),prime(s[i](1)),1.0);
        }
        CurrMPO1.setA((i+1),MPOTensor1);
        CurrMPO2.setA((i+1),MPOTensor2);
    }
    
    // call resp. position and apply operator
    rho.position(m+1);
    Args my_args = Args("Cutoff", 1e-30, "Maxm", 5000);
    MPS *res = new MPS;
    *res = MPS(N);
    MPS rho_temp1 = lv::apply_nonherm_MPO(CurrMPO1, rho, *res, my_args);
    delete res;
    
    MPS *res1 = new MPS;
    *res1 = MPS(N);
    MPS rho_temp2 = lv::apply_nonherm_MPO(CurrMPO2, rho, *res1, my_args);
    delete res1;
    
    // contract rho_MPS
    ITensor temp1 = rho_temp1.A(1);
    ITensor temp2 = rho_temp2.A(1);
    for (int k = 1; k<N; k++)
    {
        temp1 *= rho.A(k);
        temp1 *= rho_temp1.A(k+1);
        temp2 *= rho.A(k);
        temp2 *= rho_temp2.A(k+1);
    }
    temp1 *= rho.A(N);
    temp2 *= rho.A(N);
    
    // read rel_curr and scale it with trace and spin1/2 picture
    rel_curr = sqrt(temp1.real()) + sqrt(temp2.real());
    rel_curr *= 0.5*pow(2.0,N);

    return rel_curr; 
}

        
