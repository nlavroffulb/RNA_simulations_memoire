#pragma once
#include<vector>
#include "math_functions.h"
#include <fstream>
#include<string>
#include<iostream>

void sample_jump_direction(std::vector<double> &jump, double bond_distance);
std::vector<double> crankshaft_insertion(std::vector<double> &initial, std::vector<double> &N_position);
std::vector<double> rejection_sample(std::vector<double> &initial, std::vector<double> &N_position, int number_of_segments);
double ideal_chain_pdf(double &e2e_distance, int n_segments);
double segment_grow_prob(double new_e2e, double old_e2e, int l_minus_i);
bool overstretch(int l_minus_i, double &new_e2e);
double sample_e2e( double init_e2e);
std::vector<std::vector<double>> grow_chain(std::vector<double> &starting_end, 
	std::vector<double> &ending_end, int l,
	std::vector<std::vector<double>>& excluded_volumes);

std::vector<std::vector<double>> grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l);

std::vector<std::vector<double>> random_walk(std::vector<double> &starting_end, int segments, std::vector<std::vector<double>>& excluded_volumes);






