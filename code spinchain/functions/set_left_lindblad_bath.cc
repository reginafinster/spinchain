#include "itensor/all.h"
#include "header.h"

using namespace itensor;

void liouvillian::set_left_lindblad_bath (ITensor& l_tensor, const double& gamma_1_pl, const double& gamma_1_mi)
{
    // get indices
    IndexSet s = l_tensor.inds();
        
    // set diagonal entries
    for (int i = 2; i<=3; i++)
    {
        for (int j = 1; j<=4; j++)
        {
            l_tensor.set(s[0](i),s[1](j),s[2](i),s[3](j), -2.0*gamma_1_pl);
        }
    }
    
    for (int j = 1; j<=4; j++)
    {
        l_tensor.set(s[0](4),s[1](j),s[2](4),s[3](j), -4.0*gamma_1_pl);
        
    // set off-diagonal entries
        l_tensor.set(s[0](4),s[1](j),s[2](1),s[3](j), 4.0*gamma_1_mi);
    }
    
    return;
    
}
