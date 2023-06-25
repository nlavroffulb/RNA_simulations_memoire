#pragma once
#include "polymer_class.h"


void hairpin_sim(polymer* p0,int NMC);

void determine_g0(polymer* p0, int NMC);

void linker_weight_sim(polymer* p0, int NMC);

void varying_helix_length_sim(polymer* p0, int NMC);

void vary_dangle_sim(polymer* p0, int NMC);

void kissing_hairpin_swivel(polymer* p0, int NMC, bool swivel_allowed);



//void hairpin_sim_vary_WC_pot(int NMC, int N_monomers, bool rosenbluth);


