#include "itensor/all.h"
#include "header.h"  

namespace lv = liouvillian;
    
void liouvillian::fill_MPOs_Suz_Trot_2nd_order(MPO& U_odd, MPO& U_even, MPO& U_odd_thalf, MPO& U_even_thalf, const std::vector<Index>& s, const double& gamma_1_pl, const double& gamma_1_mi, const double& gamma_N_pl, const double& gamma_N_mi, const double& J, const double& Delta, const std::vector<double>& h, const double& dt, const int& include_bath_tensors)
{   
    // get number of sites
    int N = s.size();
    double Nmult = 1.0;
    ITensor left_bath_tensor = ITensor();
    ITensor right_bath_tensor = ITensor();
    ITensor l_tensor = ITensor();
    // first and last Tensor of the MPO contain the baths. Else the MPO is filled 
    // only with Hamiltonian parts
    for (int i=0; i<(N-1); i++)
    {
        // if we look at two sites only, we need to include both bath-tensors in one step
        if (N==2)
        {
           // set entries for Hamiltonian:
            l_tensor = ITensor(s[i], s[i+1], prime(s[i]), prime(s[i+1]));
            lv::set_hamilton(l_tensor, h[i], h[i+1], J, Delta);
            // factor 0.25 for spin-half-picture
            //l_tensor = 0.25*l_tensor;
            // set entries for both Lindblad-baths if switch for it is true
            if (include_bath_tensors) 
            {
                left_bath_tensor = ITensor(s[i], s[i+1], prime(s[i]), prime(s[i+1]));
                lv::set_left_lindblad_bath(left_bath_tensor, gamma_1_pl, gamma_1_mi);
                left_bath_tensor = Nmult*left_bath_tensor;
                
                right_bath_tensor = ITensor(s[i], s[i+1], prime(s[i]), prime(s[i+1]));
                lv::set_right_lindblad_bath(right_bath_tensor, gamma_N_pl, gamma_N_mi);
                right_bath_tensor = Nmult*right_bath_tensor;
                l_tensor = l_tensor + left_bath_tensor + right_bath_tensor;
            }
            
        }
        //loop through chain otherwise
        else
        {
            if (i==0)
            {
                // set entries for Hamiltonian:
                l_tensor = ITensor(s[i], s[i+1], prime(s[i]), prime(s[i+1]));
                lv::set_hamilton(l_tensor, h[i], h[i+1], J, Delta);
                // factor 0.25 for spin-half-picture
                //l_tensor = 0.25*l_tensor;
                // set entries for left Lindblad-bath if switch for it is true
                if (include_bath_tensors) 
                {
                    left_bath_tensor = ITensor(s[i], s[i+1], prime(s[i]), prime(s[i+1]));
                    lv::set_left_lindblad_bath(left_bath_tensor, gamma_1_pl, gamma_1_mi);
                    //left_bath_tensor = pow(2,(N*0.5))*left_bath_tensor;
                    left_bath_tensor = Nmult*left_bath_tensor;
                    l_tensor = l_tensor + left_bath_tensor;
                }
            }
                    
            else if (i==(N-2))
            {
                // set entries for Hamiltonian:
                l_tensor = ITensor(s[i], s[i+1],  prime(s[i]), prime(s[i+1]));
                lv::set_hamilton(l_tensor, h[i], h[i+1], J, Delta);
                // factor 0.25 for spin-half-picture
                //l_tensor = 0.25*l_tensor;
                // set entries for right Lindblad-bath if switch for it is true
                if (include_bath_tensors) 
                {
                    right_bath_tensor = ITensor(s[i], s[i+1],  prime(s[i]), prime(s[i+1]));
                    lv::set_right_lindblad_bath(right_bath_tensor, gamma_N_pl, gamma_N_mi);
                    //right_bath_tensor = pow(2,(N*0.5))*right_bath_tensor;
                    right_bath_tensor = Nmult*right_bath_tensor;
                    l_tensor = l_tensor + right_bath_tensor;
                }
            }
            else
            {
                l_tensor = ITensor(s[i], s[i+1], prime(s[i]), prime(s[i+1]));
                lv::set_hamilton(l_tensor, h[i], h[i+1], J, Delta);
                // factor 0.25 for spin-half-picture
                //l_tensor = 0.25*l_tensor;
            }
        }
        
        // write U with Taylor in second order for approximation of e-function
        // create identity tensors in both subspaces for Taylor approximation
        ITensor identity_temp_1 = ITensor(s[i], prime(s[i]));
        ITensor identity_temp_2 = ITensor(s[i+1], prime(s[i+1]));
        lv::set_identity(identity_temp_1);
        lv::set_identity(identity_temp_2);
        
        //copy tensor for step below
        ITensor l_tensor_temp = l_tensor;
        
        // make decomposition with full timestep dt
        l_tensor = l_tensor*dt;
        
        // Taylor 
        ITensor L2  = mapprime((l_tensor * prime(l_tensor)),2,1);
        ITensor L3  = mapprime((L2 * prime(l_tensor)),2,1);
        
        ITensor U_temp = identity_temp_1 * identity_temp_2 + l_tensor + (0.5 * L2) + (L3 / 6.0);
        
        // apply SVD to produce two tensors, each acting just on one single site:
        ITensor L = ITensor(s[i], prime(s[i])); // acts on site i
        ITensor R = ITensor(s[i+1], prime(s[i+1])); // acts on site i+1
        ITensor S;
        svd(U_temp,L,S,R);
        L = L * S;
        
        // write both tensors into even/odd time-evolution MPO
        if((i+1)%2==1)
        {
            U_odd.setA(i+1,L);
            U_odd.setA(i+2,R);
        }
        else
        {
            U_even.setA(i+1,L);
            U_even.setA(i+2,R);
        }
        
        //make decomposition with half timestep dt/2
        l_tensor = l_tensor_temp;
        l_tensor = l_tensor*dt*0.5;
        
        // Taylor 
        L2  = mapprime((l_tensor * prime(l_tensor)),2,1);
        L3  = mapprime((L2 * prime(l_tensor)),2,1);
        
        ITensor U_temp_thalf = identity_temp_1 * identity_temp_2 + l_tensor + (0.5 * L2) + (L3 / 6.0);
        
        // apply SVD to produce two tensors, each acting just on one single site:
        ITensor L_thalf = ITensor(s[i], prime(s[i])); // acts on site i
        ITensor R_thalf = ITensor(s[i+1], prime(s[i+1])); // acts on site i+1
        ITensor S_thalf;
        svd(U_temp_thalf,L_thalf,S_thalf,R_thalf);
        L_thalf = L_thalf * S_thalf;
        
        // write both tensors into even/odd time-evolution MPO
        if((i+1)%2==1)
        {
            U_odd_thalf.setA(i+1,L_thalf);
            U_odd_thalf.setA(i+2,R_thalf);
        }
        else
        {
            U_even_thalf.setA(i+1,L_thalf);
            U_even_thalf.setA(i+2,R_thalf);
        }
        
    }// EOF loop over sites
        
    // fill remaining entries of odd and even time-evolution operators
    // for full timestep evolution
    if(N%2==0)
    {
        ITensor identity_1 = ITensor(s[0], prime(s[0]));
        lv::set_identity(identity_1);
        U_even.setA(1, identity_1);
        ITensor identity_N = ITensor(s[N-1], prime(s[N-1]));
        lv::set_identity(identity_N);
        U_even.setA(N, identity_N);
    }
    else
    {
        ITensor identity_1 = ITensor(s[0], prime(s[0]));
        lv::set_identity(identity_1);
        U_even.setA(1, identity_1);
        ITensor identity_N = ITensor(s[N-1], prime(s[N-1]));
        lv::set_identity(identity_N);
        U_odd.setA(N, identity_N);
    }
    U_even.orthogonalize();
    U_odd.orthogonalize();
    
    if(N%2==0)
    {
        ITensor identity_1 = ITensor(s[0], prime(s[0]));
        lv::set_identity(identity_1);
        U_even_thalf.setA(1, identity_1);
        ITensor identity_N = ITensor(s[N-1], prime(s[N-1]));
        lv::set_identity(identity_N);
        U_even_thalf.setA(N, identity_N);
    }
    else
    {
        ITensor identity_1 = ITensor(s[0], prime(s[0]));
        lv::set_identity(identity_1);
        U_even_thalf.setA(1, identity_1);
        ITensor identity_N = ITensor(s[N-1], prime(s[N-1]));
        lv::set_identity(identity_N);
        U_odd_thalf.setA(N, identity_N);
    }
    // for full timestep evolution
    U_even_thalf.orthogonalize();
    U_odd_thalf.orthogonalize();
    
    return;
}
