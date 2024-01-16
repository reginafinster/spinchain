#include "itensor/all.h"
#include "header.h"

using namespace itensor;

void liouvillian::set_lindblad(ITensor& l_tensor, const double& gamma_pl, const double& gamma_mi)
{
    // get indices
    IndexSet s = l_tensor.inds();
        
    // set entries
    l_tensor.set(s[0](2),s[1](2), -2.0*gamma_pl);
    l_tensor.set(s[0](3),s[1](3), -2.0*gamma_pl);
    l_tensor.set(s[0](4),s[1](4), -4.0*gamma_pl);
    l_tensor.set(s[0](4),s[1](1), +4.0*gamma_mi);
    
    return;
    
}
