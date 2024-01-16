#include "itensor/all.h"
#include "header.h"

using namespace itensor;

void liouvillian::set_right_lindblad_bath (ITensor& l_tensor, const double& gamma_N_pl, const double& gamma_N_mi)
{
    // get indices
    IndexSet s = l_tensor.inds();
        
    for (int i = 1; i<=4; i++)
    {
        // set diagonal entries
        for (int j = 2; j<=3; j++)
        {
            l_tensor.set(s[0](i),s[1](j),s[2](i),s[3](j), -2.0*gamma_N_pl);
        }
        
        l_tensor.set(s[0](i),s[1](4),s[2](i),s[3](4), -4.0*gamma_N_pl);
        
        // set off-diagonal entries
        l_tensor.set(s[0](i),s[1](4),s[2](i),s[3](1), 4.0*gamma_N_mi);
    }

    
    return;
    
}
