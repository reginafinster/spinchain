#include "itensor/all.h"
#include "header.h"  

namespace lv = liouvillian;

double liouvillian::calc_abs_current_new(MPS& rho_MPS, const std::vector<Index>& s)
{
    // get size of MPS
    int N = rho_MPS.N();
    
    double abs_current = 0.0;
    double rel_current = 0.0;
    
    for (int i=0; i<(N-1); i++)
    {
        rel_current = lv::calc_rel_current_new(rho_MPS, i, s);
        abs_current += rel_current;
    }
    return abs_current;

}
