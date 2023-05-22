#include "helix_struct.h"

helix_struct::helix_struct(int o,std::vector<double> &u_vec, std::vector<double> &v_vec, std::vector<int> &structure, int growth_direction) {
	origin = o;
	u = u_vec;
	v = v_vec;
	monomer_indices = structure;
	alpha = growth_direction;
	//std::cout << std::endl;
}

helix_struct::helix_struct(std::vector<int> monomers, std::vector<double>& v_vec, std::vector < std::vector<double>> r_centres, int a, int b) {
	monomer_indices = monomers;
	v = v_vec;
	running_centres = r_centres; 
	alpha = a; 
	beta = b;
	set_growth_direction();
	std::vector<std::vector<double>> rcs(2);
	running_centres = rcs;

}
helix_struct::helix_struct(std::vector<int>& monomers, std::vector<double>& v_vec, int a, int b) {
	monomer_indices = monomers;
	v = v_vec;
	alpha = a;
	beta = b;
	set_growth_direction();
	std::vector<std::vector<double>> rcs(2);
	running_centres = rcs;

}
void helix_struct::set_growth_direction() {
	gamma = abs((1 - alpha) - (1 - beta));
}
void helix_struct::set_running_x_strand_centre(std::vector<double> rxc) {
	running_centres[0] = rxc;
}
void helix_struct::set_running_y_strand_centre(std::vector<double> ryc) {
	running_centres[1] = ryc;
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
void helix_struct::change_u(std::vector<double> u_vec) {
	u = u_vec;
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
int helix_struct::get_growth_direction() {
	return gamma;
}
void helix_struct::set_extension(std::vector<int> big_struct) {
	extension = big_struct;
}
void helix_struct::set_extendable(bool ex) {
	extendible = ex;
}

void helix_struct::set_origin(int o) {
	origin = o;
}
void helix_struct::extend(int o, std::vector<double>& u_vec, std::vector<double> &v_vec) {
	origin = o;
	u = u_vec;
	v = v_vec;

	monomer_indices = extension;
	//extendable = false;
	//large = true;
	//extension.clear();
}

void helix_struct::extend(std::vector<int> &x10sion) {
	monomer_indices = x10sion;
}

void helix_struct::shorten(int o, std::vector<double> &u_vec, std::vector<double> &v_vec) {
	origin = o;
	u = u_vec;
	v = v_vec;
	//extendable = true;
	//large = false;

}
void helix_struct::shorten(std::vector<int>& x10sion) {
	monomer_indices = x10sion;
}

void helix_struct::change_u_v(std::vector<double> &new_u, std::vector<double>& new_v) {
	u = new_u, v = new_v;
}

bool helix_struct::extended() {
	if (abs(monomer_indices[1] - monomer_indices[0]) + 1 ==standard_struct_size) {
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
int helix_struct::get_last_y() {
	int alpha{ get_alpha() }, beta{ get_beta() }, k{static_cast<int>(pow(-1,alpha))};

	return monomer_indices[3 * alpha + k * (beta + 2)];
}

int helix_struct::get_mon_before_origin() {
	int c{ 0 };
	int kappa{ static_cast<int>(pow(-1,get_alpha() + get_beta())) };

	int target{ origin - kappa };
	std::vector<int> reaction_region{ monomer_indices };
	int size{ abs(reaction_region[1] - reaction_region[0]) + 1 };
	std::vector<int> temp(2 * size);

	for (int i{ 0 }; i < size; i++) {
		temp[i] = reaction_region[0] + i;
	}
	for (int i{ 0 }; i < size; i++) {
		temp[i + size] = reaction_region[2] + i;
	}
	reaction_region = temp;


	for (auto i : reaction_region) {
		if (i == target) {
			return -1;
		}
	}
	return target;

}


std::vector<double> helix_struct::get_rc(int index) {
	return running_centres[index];
}
void helix_struct::set_rc(int index, std::vector<double> new_rc) {
	running_centres[index] = new_rc;

}