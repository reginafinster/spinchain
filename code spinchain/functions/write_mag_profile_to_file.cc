#include "itensor/all.h"
#include "header.h"  

namespace lv = liouvillian;

void liouvillian::write_mag_profile_to_file(ITensor& rho, const double& corr_factor, const double& timestep, FILE& file_for_mag_profile)
{
    // get indices and size of rho
    IndexSet s = rho.inds();
    int N = s.r();
    
    for (int i=0; i<N; i++)
    {
        double site_mag = lv::calc_single_site_mag(rho, i);
        fprintf(file_for_mag_profile, "%.5f \t\t %.i \t\t %.5f \n", timestep, i, (site_mag/corr_factor));
    }

}
