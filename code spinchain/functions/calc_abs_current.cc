#include "itensor/all.h"
#include "header.h"  

namespace lv = liouvillian;

double liouvillian::calc_abs_current(ITensor& rho, const double& J)
{
    // get indices and size of rho
    IndexSet s = rho.inds();
    int N = s.r();
    
    double abs_current = 0.0;
    double rel_current = 0.0;
    
    for (int i=1; i<N; i++)
    {
        rel_current = lv::calc_rel_current(rho, i, J);
        abs_current += rel_current;
    }
    return abs_current;

}
