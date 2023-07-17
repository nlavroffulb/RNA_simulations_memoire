// the last dance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include "math_functions.h"
#include "helix_generation.h"
#include "simulations.h"
#include <fstream>
#include <chrono>


int main()
{
    srand(time(NULL));
    std::int64_t init{ rand()};
    init = -1*init;
    ran2(&init);

    int NMC{ 50000 }, i{ 0 }, N_monomers(30), NMC2{ 10000 };
    std::cout << "NMC is " << NMC << std::endl;
    //polymer p0(N_monomers, true);


    //p0.set_link_acceptance_prefactor(1.0);
    // reference region to calculate prefactor
    std::vector<int> regions{ 9,11,N_monomers - 3,N_monomers - 1 };


    // test the convergence when the WCA potential is 0.
    // then the only factor should be the ideal chain configuration. that's the configurational cost. should see rapid convergence.
    // change length of linker or what? 


    //////////////////////////////////////////////////////////////////////////////////////////
    // determine appropriate weight factor
    //p0.set_link_acceptance_prefactor(0.01);
    //p0.set_search_results(regions);
    //determine_g0(&p0, NMC, 2);
    //correlation_length(p0.get_link_acceptance_prefactor(), N_monomers, regions, NMC);
    ////// varying dangle
    //double acc_prefac{ p0.get_link_acceptance_prefactor() };
    //std::cout << "varying dangle length" << std::endl;
    //vary_dangle_sim(acc_prefac, N_monomers, NMC2, regions,"Results/varying dangle sims/");
    //////////////////////////////////////////////////////////////////////////////////////
    //std::vector<int> regions{ 0,2,N_monomers - 3,N_monomers - 2 };
    //regions = { 0,2,N_monomers - 3,N_monomers - 1 };
    //polymer p1(N_monomers, true);
    //p1.set_link_acceptance_prefactor(0.01);
    //p1.set_search_results(regions);
    //determine_g0(&p1, NMC, 100);
    //correlation_length(p1.get_link_acceptance_prefactor(), N_monomers, regions, NMC);
    //vary_WCA_potential(p1.get_link_acceptance_prefactor(), N_monomers, NMC2, regions, "Results/WCA_potential/");
    ////////////////////////////////////////////////////////////////////////////////////////
    //std::vector<std::vector<int>> structures{ { 2,4,23,25 }, { 15,17,49,51 }, { 40,42,57,59 } };
    //int chain_length{ 60 };

    //polymer ref_p(chain_length, true);
    //ref_p.set_link_acceptance_prefactor(0.01);
    //p0.set_search_results(structures[1]);
    //determine_prefactor(&ref_p, NMC, 100);
    //correlation_length(ref_p.get_link_acceptance_prefactor(), chain_length, structures[1], NMC);
    //double reference_prefactor{ ref_p.get_link_acceptance_prefactor() };

    //polymer KH_OFF(chain_length, true);
    //KH_OFF.set_link_acceptance_prefactor(reference_prefactor);
    //KH_OFF.set_search_results(structures);
    ////determine_g0(&KH, NMC,100);
    ////std::cout << "past g0" << std::endl;
    //kissing_hairpin_swivel(&KH_OFF, NMC, false);

    //polymer KH_ON(chain_length, true);
    //KH_ON.set_search_results(structures);
    //KH_ON.set_link_acceptance_prefactor(reference_prefactor);
    //std::cout << "past KH with swivel off" << std::endl;
    //kissing_hairpin_swivel(&KH_ON, NMC, true);
    //////////////////////////////////////////////////////////////////////////////////////////
    //varying_helix_length(N_monomers, NMC, "Results/varying_helix_length_sims/");
    //////////////////////////////////////////////////////////////////////////////////////////
    //vary_dangle_sim(N_monomers, NMC, "Results/varying_dangle_sims/");
    testing_convergence(NMC, "Results/");

    //polymer p0(N_monomers, true);
    //std::vector<int> dh{ 5,7,27,29 };

    //int alpha, beta, struct_index;
    //std::vector<int> s;
    //std::vector<std::vector<int>> growth_limits;
    //int s_index;
    //for (int i = 0; i < 10; i++)
    //{

    //    p0.sample_link_region(s, alpha, beta, struct_index);
    //    std::cout << "alpha beta " << alpha << beta << std::endl;
    //    p0.link_growth_limits(dh, alpha, beta, growth_limits);
    //    std::cout << "link limits" << std::endl;
    //    print_2d_int_vec(growth_limits);
    //    std::cout << "unlink limits" << std::endl;
    //    p0.unlink_growth_limits(dh, alpha, beta, growth_limits);
    //    print_2d_int_vec(growth_limits);

    //    p0.link(s, alpha, beta, struct_index);
    //    p0.link_update(struct_index);
    //    p0.neighbouring_linkers();
    //    p0.reset_positions();



    //}


    //// get data
    //determine_g0(&p0, NMC,100);
    //////p0.set_link_acceptance_prefactor(0, 0.01);
    //
    //std::cout << "varying helix length" << std::endl;
    //p0.set_search_results({ 0,5,N_monomers - 6,N_monomers-1 });
    //varying_helix_length_sim(&p0, NMC2);



    //double acc_prefac{ p0.get_link_acceptance_prefactor() };
    //vary_dangle_sim(acc_prefac, N_monomers, NMC2);

    //p0.set_search_results({ 0,2,N_monomers - 3,N_monomers - 1 });
    //determine_g0(&p0, NMC,100);
    //old_new_weights_sim(&p0, NMC2);


    //polymer KH(60, true);
    //KH.set_link_acceptance_prefactor(50);

    //std::cout << "starting new simulation" << std::endl;
    //std::vector<std::vector<int>> structures{ { 2,4,23,25 }, { 15,17,49,51 }, { 40,42,57,59 } };
    //KH.set_search_results(structures);
    ////determine_g0(&KH, NMC,100);
    ////std::cout << "past g0" << std::endl;
    //kissing_hairpin_swivel(&KH, NMC, false);
    //std::cout << "past KH with swivel off" << std::endl;
    //kissing_hairpin_swivel(&KH, NMC, true);






}
