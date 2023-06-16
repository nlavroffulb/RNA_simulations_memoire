#pragma once

// we will write different versions of metropolis hastings for different situations

void link_unlink_hairpin_sim(int NMC, int N_monomers, bool rosenbluth);

void hairpin_sim_varying_length(int NMC, int N_monomers, bool rosenbluth);

void hairpin_sim_vary_WC_pot(int NMC, int N_monomers, bool rosenbluth);

void covid_simulation(int NMC);

void kissing_hairpin_simulation(int NMC);

