#pragma once
#include <cmath>
#include <string>
#include<vector>
#include "helix_struct.h"


bool base_pairing_rules(char b1, char b2);
std::vector<int> position_of_region(int i, int r);
bool regions_compatible(std::string r1, std::string r2);


void sample_helix_vectors(std::vector<double>& u, std::vector<double>& v);

double helix_separation();

//re-generating polymer with structure
std::vector<double> rotate_helix(std::vector<double> &u, std::vector<double> &v, double the);
void generate(int region_length,
    std::vector<double> &v, std::vector<double> &u, std::vector<double> &x,
    std::vector<std::vector<double>>& s_x, std::vector<std::vector<double>>& s_y);

void generate(int region_length,
    std::vector<double>& v, std::vector<double>& u, std::vector<double>& x,
    std::vector<std::vector<double>>& s_x, std::vector<std::vector<double>>& s_y, std::vector<std::vector<double>> &running_cs);

void reverse_generate(int region_length, std::vector<double>& v,
    std::vector<double>& u, std::vector<double>& x, std::vector<double> y,
    std::vector<std::vector<double>>& s_x, std::vector<std::vector<double>>& s_y);

void add_to_helix(std::vector<double> &u, std::vector<double> &v, std::vector<double> &o_x, std::vector<double> end_y, int side, int extra_bp, std::vector<std::vector<double>>& s_x, std::vector<std::vector<double>> &s_y);

double side_length(int n);

double sideways_length(int n);

double consecutive_mon_sep();
double max_ring_sep();
std::string random_base();
