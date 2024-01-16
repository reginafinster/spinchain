#include "itensor/all.h"
#include "header.h"  

namespace lv = liouvillian;
    
void liouvillian::make_lindblad_MPO(MPO& Lindblad_MPO, const std::vector<Index>& s, const double& gamma_1_pl, const double& gamma_1_mi, const double& gamma_N_pl, const double& gamma_N_mi, const double& dt)
{
    int N = s.size();
    ITensor l_tensor = ITensor();
    
    for (int i=0; i<N; i++)
    {
        l_tensor = ITensor(s[i], prime(s[i]));
        if (i==0) {lv::set_lindblad(l_tensor, gamma_1_pl, gamma_1_mi);}
        else if (i==(N-1)) {lv::set_lindblad(l_tensor, gamma_N_pl, gamma_N_mi);}
        else {lv::set_identity(l_tensor);}
        
        // write U with Taylor in second order for approximation of e-function
        // create identity tensors in both subspaces for Taylor approximation
        ITensor identity_temp = ITensor(s[i], prime(s[i]));
        lv::set_identity(identity_temp);
        
        l_tensor = l_tensor*dt;
        ITensor U_temp = identity_temp + l_tensor + (0.5 * mapprime((l_tensor * prime(l_tensor)),2,1));
        
        Lindblad_MPO.setA(i+1, U_temp);
    }
    PrintData(Lindblad_MPO);
    Lindblad_MPO.orthogonalize();
    
    return;   
    
    // TODO: weiter, das hier muss ich auf eine !!! site pro Tensor reduzieren.
}
