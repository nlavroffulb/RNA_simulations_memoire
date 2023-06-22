#pragma once
#include "polymer_class.h"


void hairpin_sim(int NMC, int N_monomers, bool rosenbluth);

void determine_g0(polymer* p0, int NMC);

void linker_weight_sim(polymer* p0, int NMC);

void varying_helix_length_sim(polymer* p0, int NMC);




//void hairpin_sim_vary_WC_pot(int NMC, int N_monomers, bool rosenbluth);


