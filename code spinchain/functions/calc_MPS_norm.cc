#include "itensor/all.h"
#include "header.h"

double liouvillian::calc_MPS_norm(MPS& psi)
{
    double norm = 0.0;
    Complex norm_cplx = 0.0;
    int N_sites = psi.N();
    
    for (int i=1; i<=N_sites; i++)
    {
        psi.position(i);
        norm_cplx += (psi.A(i) * dag(psi.A(i))).cplx();
    }
    
    norm = (norm_cplx.real()/N_sites);
    
    return norm;
}
