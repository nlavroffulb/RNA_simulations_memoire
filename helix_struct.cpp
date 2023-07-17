#include "helix_struct.h"


helix_struct::helix_struct(std::vector<int> monomers, std::vector<double>& v_vec, std::vector < std::vector<double>> r_centres, int a, int b) {
	monomer_indices = monomers;
	v = v_vec;
	running_centres = r_centres; 
	alpha = a; 
	beta = b;
	running_centres = r_centres;

}
std::vector<double> helix_struct::get_u() {
	return u;
}

std::vector<double> helix_struct::get_v() {
	return v;
}

void helix_struct::change_v(std::vector<double>v_vec) {
	v = v_vec;
}

int helix_struct::get_origin() {
	return origin;
}
std::vector<int> helix_struct::get_monomers() {
	return monomer_indices;
}
std::vector<int> helix_struct::get_extension() {
	return extension;
}
void helix_struct::set_extension(std::vector<int> big_struct) {
	extension = big_struct;
}

void helix_struct::set_origin(int o) {
	origin = o;
}

void helix_struct::extend(std::vector<int> &x10sion) {
	monomer_indices = x10sion;
}

void helix_struct::shorten(std::vector<int>& x10sion) {
	monomer_indices = x10sion;
}



bool helix_struct::extended() {
	if (std::abs(monomer_indices[1] - monomer_indices[0]) + 1 ==standard_struct_size) {
		return false;
	}
	else {
		return true;
	}
	//return large;
}

bool helix_struct::fully_extended() {
	if (monomer_indices == extension) {
		return true;
	}
	else return false;
}

bool helix_struct::extendable() {
	int c{ 0 };

	for (int i{ 0 }; i < 4; i++) {
		if (monomer_indices[i] != extension[i]) {
			c++;
		}
	}
	if (c == 0) {
		return true;
	}
	else return false;
	//for(int i{0};monomer_indc)
	//if (monomer_indices == extension) {
	//	return false;
	//}
	//else return true;
}

int helix_struct::get_alpha() {
	//if (origin == monomer_indices[0] || origin == monomer_indices[1]) {
	//	return 0;
	//}
	//else {
	//	return 1;
	//}
	return alpha;

}
int helix_struct::get_beta() {
	//if (origin == monomer_indices[0] || origin == monomer_indices[3]) {
	//	return 0;
	//}
	//else {
	//	return 1;
	//}
	return beta;

}



std::vector<double> helix_struct::get_rc(int index) {
	return running_centres[index];
}
void helix_struct::set_rc(int index, std::vector<double> new_rc) {
	running_centres[index] = new_rc;

}