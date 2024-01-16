#include "itensor/all.h"
#include "header.h"

namespace lv = liouvillian;

double liouvillian::calc_single_site_mag_new(MPS& rho, const int& m, const std::vector<Index>& s)
{
    // get size of rho
    int N = rho.N();
    
    // initialize current and tensors
    ITensor temp;
    double site_mag = 0.0;
    
    // initialize
    MPO MagMPO = MPO(N);
    ITensor MPOTensor;
    
    // set up operator 
    for (int i=0; i<N; i++)
    {
        if (i == m) 
        { 
            MPOTensor = ITensor(s[i],prime(s[i])); 
            MPOTensor.set(s[i](4),prime(s[i](4)),1.0);
        }
        else 
        { 
            MPOTensor = ITensor(s[i],prime(s[i]));
            MPOTensor.set(s[i](1),prime(s[i](1)),1.0);
        }
        MagMPO.setA((i+1),MPOTensor);
    }
    
    // call resp. position and apply operator
    rho.position(m);
    Args my_args = Args("Cutoff", 1e-30, "Maxm", 5000);
    MPS *res = new MPS;
    *res = MPS(N);
    MPS rho_temp = lv::apply_nonherm_MPO(MagMPO, rho, *res, my_args);
    delete res;
    
    // contract rho_MPS
    temp = rho_temp.A(1);
    for (int k = 1; k<N; k++)
    {
        temp *= rho.A(k);
        temp *= rho_temp.A(k+1);
    }
    temp *= rho.A(N);
    
    // read magnetization and scale it with trace and spin1/2 picture
    site_mag = temp.real();
    site_mag = sqrt(site_mag);
    site_mag *= 0.5*pow(2.0,N);
    
    return site_mag; 
}

        
