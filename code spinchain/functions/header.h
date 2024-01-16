#include "itensor/all.h"

using namespace itensor;


namespace liouvillian 
{
    void set_ground_state(MPS& psi, const std::vector<Index>& s);
    void set_hamilton(ITensor& l_tensor, const double& h_one, const double& h_two, const double& J, const double& Delta);
    void set_left_lindblad_bath (ITensor& l_tensor, const double& gamma_1_pl, const double& gamma_1_mi);
    void set_right_lindblad_bath (ITensor& l_tensor, const double& gamma_N_pl, const double& gamma_N_mi);
    void set_identity(ITensor& my_identity);
    void MultiplyTwoMPO(MPO const& Aorig, MPO const& Borig, MPO& res, Args args);
    void fill_MPOs_Suz_Trot(MPO& U_odd, MPO& U_even, const std::vector<Index>& s, const double& gamma_1_plus, const double& gamma_1_minus, const double& gamma_N_plus, const double& gamma_N_minus, const double& J, const double& Delta, const std::vector<double>& h, const double& dt, const int& include_bath_tensors);
    void fill_MPOs_Suz_Trot_2nd_order(MPO& U_odd, MPO& U_even, MPO& U_odd_thalf, MPO& U_even_thalf, const std::vector<Index>& s, const double& gamma_1_plus, const double& gamma_1_minus, const double& gamma_N_plus, const double& gamma_N_minus, const double& J, const double& Delta, const std::vector<double>& h, const double& dt, const int& include_bath_tensors);
    double calc_MPS_norm(MPS& psi);
    void set_sigma_z(ITensor& sigma_z);
    double trace_check(const ITensor& rho);
    double calc_stag_mag(ITensor& rho);
    MPS apply_nonherm_MPO(MPO& L, MPS& psi, MPS& res, Args args);
    MPO multiply_two_nonherm_MPO(const MPO& MPO_A, const MPO& MPO_B);
    void set_alldown_oneup_state(MPS& psi, const std::vector<Index>& s, int k);
    double calc_single_site_mag(ITensor& rho, int i);
    void set_allup_state(MPS& psi, const std::vector<Index>& s);
    void set_steady_state(MPS& psi, const std::vector<Index>& s);
    void set_neel_state_first_up(MPS& psi, const std::vector<Index>& s);
    void set_neel_state_first_down(MPS& psi, const std::vector<Index>& s);
    void set_allup_onedown_state(MPS& psi, const std::vector<Index>& s, int k);
    bool steady_state_check(std::vector<double>& cont, const double& mag, const double& tol);
    double calc_rel_current(const ITensor& rho, const int& i, const double& J);
    void set_disorder(std::vector<double>& h, std::vector<double>& check_average, const double& disorder);
    double calc_abs_current(ITensor& rho, const double& J);
    double calc_single_site_mag_new(MPS& rho, const int& m, const std::vector<Index>& s);
    double trace_check_new(MPS& rho, const std::vector<Index>& s);
    double calc_abs_current_new(MPS& rho_MPS, const std::vector<Index>& s);
    double calc_rel_current_new(MPS& rho, const int& m, const std::vector<Index>& s);
}
