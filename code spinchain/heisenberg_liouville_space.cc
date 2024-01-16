#include <stdio.h>
#include <time.h>
#include <complex>
#include <math.h>
#include <sys/stat.h> 
#include <sys/types.h> 
#include "itensor/all.h"
#include "functions/header.h"


using namespace itensor;
namespace lv = liouvillian;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Heisenberg XXZ model, spin chain with N sites, Lindblad-bath for 1st and Nth site. Basis: Pauli-Matrices     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//this program calculates the time evolution of the density operator solving the von-Neumann-equation
//with a Liouville-operator for a chain of N spin-1/2-sites with a next-neighbor-interaction (Heisenberg-Chain)
//and Lindblad-Baths at both ends of the chain.
//The basis used are the Pauli-matrices. 



int main(int argc, char* argv[]) 
{
    if(argc!=2) 
    {
        printfln("Specify input file");
        return 0;
    }
    
    // set up parameters from input file
    // params for unitary evolution
    auto input = InputGroup(argv[1],"input");
    int N = input.getInt("N", 4); 
    double J = input.getReal("J", 1.0); 
    double Delta = input.getReal("Delta", 1.0);

    double my_cutoff = input.getReal("cutoff", 1E-14); // cutoff of Schmidt values
    int my_maxm = input.getInt("maxm", 5000); // maximal number of Schmidt values 
    
    // params for Lindblad baths
    int include_bath_tensors = input.getInt("include_bath_tensors", 0);
    double mu_lower = input.getReal("mu", 0.5); // thermodynamic potential
    double gamma_1_plus = input.getReal("gamma", 1.0); // inscattering on site 1
    double gamma_1_minus = input.getReal("gamma", 1.0); // outscattering on site 1
    double gamma_N_plus = input.getReal("gamma", 1.0); // inscattering on site N
    double gamma_N_minus = input.getReal("gamma", 1.0); // outscattering on site N
    
    // param for disorder
    double disorder = input.getReal("disorder", 0.0);
    
    // params for time evolution
    int j = input.getInt("j", 1); // output is written every jth timestep
    double dt = input.getReal("dt", 0.1); // timestep
    double T = input.getReal("T", 0.2); // total time 
    
    // params for initial state: only one must be true, all other false
    int ground_state = input.getInt("ground_state", 0);
    int allup_state = input.getInt("allup_state", 0);
    int steady_state = input.getInt("quasi_steady_state", 0);
    int neel_state_first_up = input.getInt("neel_state_first_up", 0);
    int neel_state_first_down = input.getInt("neel_state_first_down", 0);
    int alldown_oneup = input.getInt("alldown_oneup", 0);
    int allup_onedown = input.getInt("allup_onedown", 0);
    int site_for_up_or_down = input.getInt("site_for_up_or_down", 0); // must be 0,1,2,3
    
    // params for calculating of magnetization: only one must be true, all other false
    int calc_mag_profile = input.getInt("calc_mag_profile", 1);
    int calc_current_profile = input.getInt("calc_current_profile", 1);
    int calc_entanglement_profile = input.getInt("calc_entanglement_profile", 1);
    
    // params for calc of steady state and output
    int check_steady_state = input.getInt("check_steady_state", 0);
    int trotter_2nd_order = input.getInt("trotter_2nd_order", 0);
    int trotter_3rd_order = input.getInt("trotter_3rd_order", 0);
    unsigned int min_const_timesteps_for_steady_state = input.getInt("min_const_timesteps_for_steady_state", 300);
    double tol = input.getReal("tolerance_for_steady_state", 0.001);
    
    int use_contraction_calc = input.getInt("use_contraction_calc", 0);
    int use_operator_calc = input.getInt("use_operator_calc", 0);
    int timesteps_before_calc = input.getInt("timesteps_before_calc", 0);
    
    // params for average
    int no_average = input.getInt("no_average", 1);
    
    
    // calculate additional parameters from input
    int timesteps = int(T/dt+(1e-9*(T/dt)));
    //int timesteps = int(floor(T/dt));
    
    // correction factor (because MPS norm is initialized to 1) for magnetization and current 
    double corr_factor = pow(2,(N*0.5));
    
    printf("Setting up directories and files\n\n");
    std::ostringstream dir_name;
    dir_name << "./results";
    // dir_name << "./MPS_N" << N << "_Tr2" << trotter_2nd_order << "_Tr3" << trotter_3rd_order << "_J" << J << "_G" << gamma_1_plus << "_T" << T << "_dt" << dt << "_mu" << mu_lower << "_c" << my_cutoff << "_maxm" << my_maxm << "_h" << disorder << "_av" << no_average << "_results";
    mkdir(dir_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::vector<Index> s = std::vector<Index>(N);
    for (int i=0; i<N; i++)
    {
        s.at(i) = Index(nameint("s",i), 4, Site);
    }
        
    std::vector<double> check_average = std::vector<double>(N);
    // für die Mittelung: pro Zeitschritt ein Eintrag, um über verschiedene Realisierungen zu mitteln.
    double av_stag_mag_trans [timesteps+1] = {};
    double av_abs_curr_trans [timesteps+1] = {};
    
    // Mittelung für steady state Werte
    double av_rel_cur = 0.0;
    double av_abs_cur = 0.0;
    double av_maxbond_size = 0.0;
    double av_avebond_size = 0.0;
    double av_mag = 0.0;
    double av_norm = 0.0;
    double av_trace = 0.0;
    bool check = 0;
    
    // array für die Berechnung der Standardabweichung
    double curr_container [no_average] = {};

    for (int m=1; m<=no_average; m++)
    {
        printf("Starting average number %i of %i \n", m, no_average);
        
        // start measuring runtime for setting up mps and mpo
        time_t starttime_mpo = time(NULL);
        
        // set up rho als MPS, write initial state and orthogolanize
        MPS *rho_MPS = new MPS;
        *rho_MPS = MPS(N);

        if (ground_state) {lv::set_ground_state(*rho_MPS, s);}
        else if (allup_state) {lv::set_allup_state(*rho_MPS, s);}
        else if (steady_state) {lv::set_steady_state(*rho_MPS, s);}
        else if (neel_state_first_up) {lv::set_neel_state_first_up(*rho_MPS, s);}
        else if (neel_state_first_down) {lv::set_neel_state_first_down(*rho_MPS, s);}
        else if (alldown_oneup) {lv::set_alldown_oneup_state(*rho_MPS, s, site_for_up_or_down);}
        else if (allup_onedown) {lv::set_allup_onedown_state(*rho_MPS, s, site_for_up_or_down);}

        printf("Initializing rho_MPS in the following state: \n");
        PrintData(*rho_MPS);
        (*rho_MPS).orthogonalize();

        // initialize container for disorder:
        std::vector<double> h = std::vector<double>(N);
        // set disorder with random generator
        lv::set_disorder(h, check_average, disorder);
        
        //create subfolder
        std::ostringstream subdir_name;
        subdir_name << dir_name.str().c_str() << "/transient_results";
        mkdir(subdir_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        // set up output stream
        std::ostringstream filename_magprofile;
        filename_magprofile << subdir_name.str().c_str() << "/mag_profile" << N << ".dat";
        FILE* file_for_mag_profile;
        file_for_mag_profile = fopen(filename_magprofile.str().c_str(), "w");
        fprintf(file_for_mag_profile, "# time \t\t magnetization \n");
        
        std::ostringstream filename_currentprofile;
        filename_currentprofile << subdir_name.str().c_str() << "/current_profile" << N << ".dat";
        FILE* file_for_current_profile;
        file_for_current_profile = fopen(filename_currentprofile.str().c_str(), "w");
        fprintf(file_for_current_profile, "# time \t\t rel current \n");
        
        std::ostringstream filename_entanglementprofile;
        filename_entanglementprofile << subdir_name.str().c_str() << "/entanglement_profile" << N << ".dat";
        FILE* file_for_entanglement_profile;
        file_for_entanglement_profile = fopen(filename_entanglementprofile.str().c_str(), "w");
        fprintf(file_for_entanglement_profile, "# time \t\t size of bond index \n");
        
        std::ostringstream filename;
        filename << subdir_name.str().c_str() << "/MPS_N" << N << "_J" << J << "_D" << Delta << "_G" << gamma_1_plus << "_mu" << mu_lower << "_dt" << dt << "_cutoff" << my_cutoff << "_maxm" << my_maxm << ".dat";
        FILE* file;
        file = fopen(filename.str().c_str(), "w");
        fprintf(file, "# time \t\t abs current \t  magnetization \t trace \t\t norm \t\t averageBond \t maxBond \t isOrthCheck \n");
            
        // calculate inscattering and outscattering rates: 
        double L_N_plus = gamma_N_plus * (1.0 - mu_lower);
        double L_N_minus = gamma_N_minus * (1.0 + mu_lower);
        
        double L_1_plus = gamma_1_plus * (1.0 + mu_lower);
        double L_1_minus = gamma_1_minus * (1.0-mu_lower);
        
        double gamma_N_pl_eff = 0.5*(L_N_plus + L_N_minus);
        double gamma_N_mi_eff = 0.5*(L_N_plus - L_N_minus);
        
        double gamma_1_pl_eff = 0.5*(L_1_plus + L_1_minus);
        double gamma_1_mi_eff = 0.5*(L_1_plus - L_1_minus);
        
        printf("Rates for coupling to the baths are: L_N_plus = %.5f, L_N_minus = %.5f, L_1_plus = %.5f, L_1_minus = %.5f \n", L_N_plus, L_N_minus, L_1_plus, L_1_minus);
        printf("gamma_N_plus_eff = %.5f, gamma_N_minus_eff = %.5f, gamma_1_plus_eff = %.5f, gamma_1_minus_eff = %.5f \n", gamma_N_pl_eff, gamma_N_mi_eff, gamma_1_pl_eff, gamma_1_mi_eff);
        
        //set up Liouvillian as MPO using the Suzukki-Trotter-Decomposition
        MPO U_MPO = MPO(N);
        Args my_args = Args("Cutoff", my_cutoff, "Maxm", my_maxm);
        
        if (trotter_3rd_order)
        {
            // additional arguments for third order decomposition
            MPO U_odd = MPO(N);
            MPO U_even = MPO(N);
            MPO U_odd_1 = MPO(N);
            MPO U_even_1 = MPO(N);
            MPO U_odd_thalf_1 = MPO(N);
            MPO U_even_thalf_1 = MPO(N);
            MPO U_MPO_1 = MPO(N);
            
            MPO U_odd_3 = MPO(N);
            MPO U_even_3 = MPO(N);
            MPO U_odd_thalf_3 = MPO(N);
            MPO U_even_thalf_3 = MPO(N);
            MPO U_MPO_3 = MPO(N);
            
            // calculate two different values for tau 
            
            double tau_1 = pow((4.0 - pow(4.0,-(1.0/3.0))),-1.0) * dt;
            double tau_3  = dt - 4*tau_1;
            
            // set up first/second part of U
            lv::fill_MPOs_Suz_Trot_2nd_order(U_odd_1, U_even_1, U_odd_thalf_1, U_even_thalf_1, s, gamma_1_pl_eff, gamma_1_mi_eff, gamma_N_pl_eff, gamma_N_mi_eff, J, Delta, h, tau_1, include_bath_tensors);
            U_MPO_1 = lv::multiply_two_nonherm_MPO(U_even_1, U_odd_thalf_1);
            U_MPO_1 = lv::multiply_two_nonherm_MPO(U_odd_thalf_1, U_MPO_1);
            
            // set up third part of U
            lv::fill_MPOs_Suz_Trot_2nd_order(U_odd_3, U_even_3, U_odd_thalf_3, U_even_thalf_3, s, gamma_1_pl_eff, gamma_1_mi_eff, gamma_N_pl_eff, gamma_N_mi_eff, J, Delta, h, tau_3, include_bath_tensors);
            U_MPO_3 = lv::multiply_two_nonherm_MPO(U_even_3, U_odd_thalf_3);
            U_MPO_3 = lv::multiply_two_nonherm_MPO(U_odd_thalf_3, U_MPO_3);
            
            // multiply all Hamiltonians to one
            U_MPO = lv::multiply_two_nonherm_MPO(U_MPO_1, U_MPO_1);
            U_MPO = lv::multiply_two_nonherm_MPO(U_MPO_3, U_MPO);
            U_MPO = lv::multiply_two_nonherm_MPO(U_MPO_1, U_MPO);
            U_MPO = lv::multiply_two_nonherm_MPO(U_MPO_1, U_MPO);
            
        }
        
        if (trotter_2nd_order)
        {
            // additional arguments for second order decomposition
            MPO U_odd = MPO(N);
            MPO U_even = MPO(N);
            MPO U_odd_thalf = MPO(N);
            MPO U_even_thalf = MPO(N);
            MPO U_MPO_temp = MPO(N);
            
            lv::fill_MPOs_Suz_Trot_2nd_order(U_odd, U_even, U_odd_thalf, U_even_thalf, s, gamma_1_pl_eff, gamma_1_mi_eff, gamma_N_pl_eff, gamma_N_mi_eff, J, Delta, h, dt, include_bath_tensors);
            
            // and contract both MPOs into a single one
            //lv::MultiplyTwoMPO(U_even, U_odd, U_MPO, my_args);
            U_MPO_temp = lv::multiply_two_nonherm_MPO(U_even, U_odd_thalf);
            U_MPO = lv::multiply_two_nonherm_MPO(U_odd_thalf, U_MPO_temp);
        }
        else 
        {
            MPO U_odd = MPO(N);
            MPO U_even = MPO(N);
            
            lv::fill_MPOs_Suz_Trot(U_odd, U_even, s, gamma_1_pl_eff, gamma_1_mi_eff, gamma_N_pl_eff, gamma_N_mi_eff, J, Delta, h, dt, include_bath_tensors);
            
            // and contract both MPOs into a single one:
            
            //lv::MultiplyTwoMPO(U_even, U_odd, U_MPO, my_args);
            U_MPO = lv::multiply_two_nonherm_MPO(U_even, U_odd);
        }
        
        time_t endtime_mpo = time(NULL);
        printf("Finished setting up MPS and MPO. Time spent here: %.1f minutes \n\n", (endtime_mpo-starttime_mpo)/60.0 );
        
        // initialize observables
        float norm = 0.0;
        float trace = 0.0;
        float mag = 0.0;
        float rel_current = 0.0;
        float abs_current = 0.0;
        int maxbond_size = 0;
        float avebond_size = 0.0;
        ITensor rho_contracted;
        
        // set up params for check for steady state:
        std::vector<double> container = std::vector<double>();
        unsigned int steady_state_counter = 0;
        int reached_steady_state = 0;
        
        // start measuring runtime for time integration
        time_t starttime_int = time(NULL); 
        
        // time loop
        for(int i = 0; i <= timesteps; i++)
        {
            if (i%10 == 0) {printf("J*dt %.2f of J*T=%.2f \n", (i*dt*J), (J*T));}
            
            if (use_contraction_calc == 1) 
            {
                // contract and output rho every jth timestep
                if(i%j == 0) 
                {
                    // check norm of rho
                    norm = lv::calc_MPS_norm(*rho_MPS);
                    check = (*rho_MPS).isOrtho();
                    
                    // contract rho
                    rho_contracted = (*rho_MPS).A(1);
                    for (int k = 2; k<=N; k++)
                    {
                        rho_contracted *= (*rho_MPS).A(k);
                    }
                    
                    trace = lv::trace_check(rho_contracted);
                    mag = lv::calc_stag_mag(rho_contracted);
                    abs_current = lv::calc_abs_current(rho_contracted, J);
                    
                    // write stag mag and abs curr into array at position of timestep and divide with number of averages only when writing to file:
                    av_stag_mag_trans[i] += mag;
                    av_abs_curr_trans[i] += abs_current;
                
                    if (calc_mag_profile) 
                    {
                        fprintf(file_for_mag_profile, "%.5f \t\t" , (i*dt*J));
                        for (int m=0; m<N; m++)
                        {
                            double site_mag = lv::calc_single_site_mag(rho_contracted, m);
                            fprintf(file_for_mag_profile, "%.5f \t\t", (site_mag/corr_factor));
                        }
                        fprintf(file_for_mag_profile, "\n");
                    }

                    if (calc_current_profile) 
                    {
                        fprintf(file_for_current_profile, "%.5f \t\t" , (i*dt*J));
                        for (int m=1; m<N; m++)
                        {
                            double site_current = lv::calc_rel_current(rho_contracted, m, J);
                            fprintf(file_for_current_profile, "%.5f \t\t", (site_current/corr_factor));
                        }
                        fprintf(file_for_current_profile, "\n");
                    }
                    
                    // write output into file
                    fprintf(file, "%.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.i \n", (i*dt*J), (av_abs_curr_trans[i]/(corr_factor*double(m))), (av_stag_mag_trans[i]/(corr_factor*double(m))), trace*pow(2.0,N), norm, check); 
                    
                    // fprintf(file, "%.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.i \n", (i*dt*J), (abs_current/corr_factor), (mag/corr_factor), trace*pow(2.0,N), norm, check); 
                    
                    // check if steady state is reached
                    if (check_steady_state)
                    {
                        if (steady_state_counter < min_const_timesteps_for_steady_state)
                        {
                            container.push_back(abs_current);
                            steady_state_counter++;
                        }
                        else 
                        {
                            container.erase(container.begin());
                            container.push_back(abs_current);
                            reached_steady_state = lv::steady_state_check(container, abs_current, tol);
                        }
                    }
                }
                
                // if steady state is reached, stop time integration
                if (reached_steady_state) {printf("reached steady state! \n"); break;}
                
                // apply liouvillian on rho
                //rho_MPS = exactApplyMPO(U_MPO, rho_MPS, my_args);
                MPS *res = new MPS;
                *res = MPS(N);
                *rho_MPS = lv::apply_nonherm_MPO(U_MPO, *rho_MPS, *res, my_args);
                delete res;
                (*rho_MPS).orthogonalize();
                
                
            } // EOF contraction loop
        
            else if (use_operator_calc == 1) 
            {
                if (i<timesteps_before_calc) 
                {
                    MPS *res = new MPS;
                    *res = MPS(N);
                    *rho_MPS = lv::apply_nonherm_MPO(U_MPO, *rho_MPS, *res, my_args);
                    delete res;
                    (*rho_MPS).orthogonalize();
                }
                
                else 
                {    
                    if (i%j == 0) // calculation of control parameters and expectation values only in this loop
                    {
                        norm = lv::calc_MPS_norm(*rho_MPS);
                        check = (*rho_MPS).isOrtho();
                        trace = (lv::trace_check_new(*rho_MPS, s))/corr_factor;
                        
                        if (calc_mag_profile) 
                        {
                            fprintf(file_for_mag_profile, "%.5f \t\t" , (i*dt*J));
                            for (int m=0; m<N; m++)
                            {
                                double site_mag = lv::calc_single_site_mag_new(*rho_MPS, m, s);
                                fprintf(file_for_mag_profile, "%.5f \t\t ", site_mag/corr_factor);
                            }
                            fprintf(file_for_mag_profile, "\n");
                        }

                        if (calc_current_profile) 
                        {
                            fprintf(file_for_current_profile, "%.5f \t\t" , (i*dt*J));
                            for (int m=0; m<(N-1); m++)
                            {
                                double site_current = lv::calc_rel_current_new(*rho_MPS, m, s);
                                fprintf(file_for_current_profile, "%.5f \t\t", (site_current/corr_factor));
                            }
                            fprintf(file_for_current_profile, "\n");
                        }
                        
                        if (calc_entanglement_profile)
                        {
                            fprintf(file_for_entanglement_profile, "%.5f \t\t" , (i*dt*J));
                            for (int m=0; m<(N-1); m++)
                            {
                                Index bond_index = commonIndex((*rho_MPS).A(m+1), (*rho_MPS).A(m+2), Link);
                                long int bond_index_size = bond_index.m();
                                fprintf(file_for_entanglement_profile, "%li \t\t", bond_index_size);
                            }
                            fprintf(file_for_entanglement_profile, "\n");
                            
                        }
                        
                        // calculate absolute current
                        // no average over transient dynamics here because of missing sign in value
                        abs_current = lv::calc_abs_current_new(*rho_MPS, s);
                        maxbond_size = maxM(*rho_MPS);
                        avebond_size = averageM(*rho_MPS);
                        
                        // write output into file
                        fprintf(file, "%.5f \t\t  %.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.i \t\t %.i \n", (i*dt*J), (abs_current/corr_factor), (mag/corr_factor), trace*pow(2.0,N), norm, avebond_size, maxbond_size, check);
                        
                        // check if steady state is reached
                        if (check_steady_state)
                        {
                            if (steady_state_counter < min_const_timesteps_for_steady_state)
                            {
                                container.push_back(abs_current);
                                steady_state_counter++;
                            }
                            else 
                            {
                                container.erase(container.begin());
                                container.push_back(abs_current);
                                reached_steady_state = lv::steady_state_check(container, abs_current, tol);
                            }
                        }
                        
                    }// EOF calculation loop
                    
                    // if steady state is reached, stop time integration
                    if (reached_steady_state) {printf("reached steady state! \n"); break;}
                    
                    MPS *res = new MPS;
                    *res = MPS(N);
                    *rho_MPS = lv::apply_nonherm_MPO(U_MPO, *rho_MPS, *res, my_args);
                    delete res;
                    (*rho_MPS).orthogonalize();
                }
                
            } // EOF calc_operator_loop
                
        }// EOF time loop
        
        fclose(file);
        fclose(file_for_mag_profile);
        fclose(file_for_current_profile);
        fclose(file_for_entanglement_profile);
    
        delete rho_MPS;
        
        // measure time
        time_t endtime_int = time(NULL);
        printf("Finished time integration. Time spent here: %.1f minutes \n\n", (endtime_int-starttime_int)/60.0 );
        
        rel_current = abs_current/double(N-1);
        
        av_maxbond_size += double(maxbond_size);
        av_avebond_size += avebond_size;
        
        av_rel_cur += rel_current;
        av_abs_cur += abs_current;
        
        av_mag += mag;
        av_norm += norm;
        av_trace += trace;
        
        curr_container[m-1] = abs_current;
        
    } // EOF average loop
    
    
    std::ostringstream filename_steadystate;
    filename_steadystate << dir_name.str().c_str() << "/steadystate.dat";
    FILE* file_steadystate;
    file_steadystate = fopen(filename_steadystate.str().c_str(), "w");

    printf("Finished loop over all realizations. \nMean applied disorder:\n" );
    for (int m=0; m<N; m++)
    {
        printf("#Average disorder at site %i: %.5f \n", (m+1), check_average[m]);
        fprintf(file_steadystate, "#Average disorder at site %i: %.5f \n", (m+1), check_average[m]);
    }
    
    //scale results 
    av_rel_cur = av_rel_cur/double(no_average);
    av_abs_cur = av_abs_cur/double(no_average);
    av_maxbond_size = av_maxbond_size/double(no_average);
    av_avebond_size = av_avebond_size/double(no_average);
    av_mag = av_mag/double(no_average);
    av_norm = av_norm/double(no_average);
    av_trace = av_trace/double(no_average);
    
    
    // calculate standard deviation
    double standard_deviation = 0.0;
    for (int m=0; m<no_average; m++)
    {
        standard_deviation += pow(curr_container[m] - av_abs_cur, 2.0);
    }
    standard_deviation /= (no_average-1);
    standard_deviation = pow(standard_deviation,0.5);
    double abs_scatter_error = standard_deviation / pow(no_average, 0.5);
    double rel_scatter_error = abs_scatter_error/fabs(av_abs_cur);
    fprintf(file_steadystate, "#Standard deviation is: %.3f. Relative error of current is: %.3f. Absolute Error is: %.3f\n", standard_deviation, rel_scatter_error, abs_scatter_error);
    
    
    // write results to file
    fprintf(file_steadystate, "# abs_current \t\t rel_current \t\t mag \t\t trace \t\t\t\t norm \t\t averageBondSize \t\t maxBondSize \t\t isOrthCheck \t\t N \t\t totaltime \t\t J \t\t Delta \t\t mu \t\t gamma \t\t cutoff \t\t maxm \t\t dt \n");
    fprintf(file_steadystate, "%.10f \t\t %.10f \t\t %.10f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.5f \t\t %.i \t\t %.i \t\t %.3f \t\t %.3f \t\t %.3f \t\t %.3f \t\t %.3f \t\t %.e \t\t %.i \t\t %.3f \n", (av_abs_cur/corr_factor), (av_rel_cur/corr_factor), (av_mag/corr_factor), av_trace*pow(2.0,N), av_norm, av_avebond_size, av_maxbond_size, check, N, T, J, Delta, mu_lower, gamma_1_plus, my_cutoff, my_maxm, dt);
    
    fclose(file_steadystate);
    return 0;
}









