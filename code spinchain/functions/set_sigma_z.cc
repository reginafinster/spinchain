#include "itensor/all.h"
#include "header.h"

void liouvillian::set_sigma_z(ITensor& sigma_z)
{
    IndexSet s = sigma_z.inds();
    
    sigma_z.set(s[0](1), s[1](4), 1.0);
    sigma_z.set(s[0](2), s[1](3), -Cplx_i);
    sigma_z.set(s[0](3), s[1](2), Cplx_i);
    sigma_z.set(s[0](4), s[1](1), 1.0);
}
