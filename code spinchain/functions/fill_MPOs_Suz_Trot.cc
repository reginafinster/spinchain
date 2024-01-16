#include "itensor/all.h"
#include "header.h"  

namespace lv = liouvillian;
    
void liouvillian::fill_MPOs_Suz_Trot(MPO& U_odd, MPO& U_even, const std::vector<Index>& s, const double& gamma_1_pl, const double& gamma_1_mi, const double& gamma_N_pl, const double& gamma_N_mi, const double& J, const double& Delta, const std::vector<double>& h, const double& dt, const int& include_bath_tensors)
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
        
        l_tensor = l_tensor*dt;
        // WARNING Im Ergebnis des Quadrats ist die Reihenfolge der Indizes vertauscht. Das stört mich. Wenn ich nichts mehr befülle, ist es kein Problem. Aber sollte ich noch was reinschreiben wollen, müsste ich das dringend beachten. 
        //ITensor U_temp = identity_temp_1*identity_temp_2 + l_tensor + (0.5 * mapprime((l_tensor * prime(l_tensor)),2,1));
        
        // Taylor in thenth order
        ITensor L2  = mapprime((l_tensor * prime(l_tensor)),2,1);
        ITensor L3  = mapprime((L2 * prime(l_tensor)),2,1);
        ITensor L4  = mapprime((L3 * prime(l_tensor)),2,1);
        ITensor L5  = mapprime((L4 * prime(l_tensor)),2,1);
        ITensor L6  = mapprime((L5 * prime(l_tensor)),2,1);
        ITensor L7  = mapprime((L6 * prime(l_tensor)),2,1);
        ITensor L8  = mapprime((L7 * prime(l_tensor)),2,1);
        ITensor L9  = mapprime((L8 * prime(l_tensor)),2,1);
        ITensor L10 = mapprime((L9 * prime(l_tensor)),2,1);
        
        ITensor U_temp = identity_temp_1 * identity_temp_2 + l_tensor + (0.5 * L2) + (L3 / 6.0) + (L4 / 24.0) + (L5 / 120.0) + (L6 / 720.0) + (L7 / 5040.0) + (L8 / 40320.0) + (L9 / 362880.0) + (L10 / 3628800.0);
        
        
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
    }// EOF loop over sites
        
    // fill remaining entries of odd and even time-evolution operators
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
    return;
}
