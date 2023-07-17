#pragma once
#include "polymer_class.h"

void hairpin_sim(polymer* p0,int NMC);

void determine_prefactor(polymer* p0, int NMC, int update_rate, std::string filename="");

double get_prefactor(std::vector<int> structure, int N_monomers, int NMC, int update_rate, std::string filename = "");

void linker_weight_sim(int N_monomers, int NMC, std::string filename = "");

void varying_helix_length(int N_monomers, int NMC, std::string filename="");


//void vary_dangle_sim(polymer* p0, int NMC, std::string filename="");
//void vary_dangle_sim(double acceptance_prefactor, int N_monomers, int NMC, const std::vector<int>& search_result, std::string filename = "");

void vary_dangle_sim(int N_monomers, int NMC, std::string filename = "");

double correlation_length(double acceptance_prefactor, int N_monomers, std::vector<int> structure, int NMC, std::string filename="");
void kissing_hairpin_swivel(polymer* p0, int NMC, bool swivel_allowed, std::string filename="");

void old_new_weights_sim(polymer* p0, int NMC, std::string filename="");

void vary_WCA_potential(double prefactor, int N_monomers, int NMC, const std::vector<int> &search_result, std::string path  = "");
//void hairpin_sim_vary_WC_pot(int NMC, int N_monomers, bool rosenbluth);

void testing_convergence(int NMC, std::string path = "");



