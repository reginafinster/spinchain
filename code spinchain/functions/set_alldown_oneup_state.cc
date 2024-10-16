#include "itensor/all.h"
#include "header.h"


void liouvillian::set_alldown_oneup_state(MPS& psi, const std::vector<Index>& s, int k)
{    
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
            // state k is up
            if (i == k)
            {
                rho_tensor = ITensor(s[i]);
                rho_tensor.set(s[i](1),0.5);
                rho_tensor.set(s[i](4),0.5);
                rho_tensor *= sqrt(2);
                psi.setA((i+1),rho_tensor);
            }
            
            // all other states down
            else
            {
                rho_tensor = ITensor(s[i]);
                rho_tensor.set(s[i](1),0.5);
                rho_tensor.set(s[i](4),-0.5);
                rho_tensor *= sqrt(2);
                psi.setA((i+1),rho_tensor);
            }
        }
    }
    return;
}
