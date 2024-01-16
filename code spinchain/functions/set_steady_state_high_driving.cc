#include "itensor/all.h"
#include "header.h"


void liouvillian::set_steady_state(MPS& psi, const std::vector<Index>& s)
{    
    // this function sets the steady state for weak driving, i.e. mu = 0.1 or mu = 0.01, and for an high inscattering left and outscattering right.
    // get size of MPS
    int N = s.size();
    int N_MPS = psi.N();
    if (N != N_MPS) 
    {
        std::cout << "mismatch in dimensions between index container and MPS" << std::endl;
    }
    
    else 
    {
        // initialize tensor
        ITensor rho_tensor = ITensor();
        
        // set tensor and fill MPS
        for (int i=0; i<N; i++)
        {
            if (i==0) // first one is up
            {
                rho_tensor = ITensor(s[i]);
                rho_tensor.set(s[i](1),0.5);
                rho_tensor.set(s[i](4),0.5);
                rho_tensor *= sqrt(2);
                psi.setA((i+1),rho_tensor);
            }
            else if (i==(N-1)) // last one is down
            {
                rho_tensor = ITensor(s[i]);
                rho_tensor.set(s[i](1),0.5);
                rho_tensor.set(s[i](4),-0.5);
                rho_tensor *= sqrt(2);
                psi.setA((i+1),rho_tensor);
            }
            else // |psi>= 1/sqrt(2)(|0>+|1>) translates to 1/2(sigma_0 + sigma_x) in Spin-Basis
            {
                rho_tensor = ITensor(s[i]);
                rho_tensor.set(s[i](1),0.5);
                rho_tensor.set(s[i](2),0.5);
                rho_tensor *= sqrt(2);
                psi.setA((i+1),rho_tensor);
            }
        }
    }
    return;
}
