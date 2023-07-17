#include <vector>
#include "math_functions.h"
#include "polymer_generation.h"
#include "helix_generation.h"
#include "helix_struct.h"
class rosenbluth_growth {
private:

	//wicks chandler parameters


	std::vector<double> weights;
	std::vector<double> energies;
	//std::vector<int> monomer_limits;


public:

	rosenbluth_growth(int length);
	rosenbluth_growth() = default;
	rosenbluth_growth(std::vector<double>& weight_vals, std::vector<double>& U_values);
	~rosenbluth_growth();
	double get_config_prob();
	double get_rosen_weight();
	double get_energy();
	double get_Ui(int index);
	double get_Wi(int index);

	std::vector<double> get_energies();
	void set_energies(std::vector<double> energy_vals);
	std::vector<double> get_weights();
	void set_weights(std::vector<double> weight_vals);

	double subsection_weight(std::vector<int> limit);
	double subsection_energy(std::vector<int> limit);
	double subsection_probability(std::vector<int> &limit);

	void modify_weights(std::vector<int> limits, std::vector<double> new_weights);
	void modify_energies(std::vector<int> limits, std::vector<double> new_energies);

	rosenbluth_growth& operator=(const rosenbluth_growth& rhs);
};

void set_WCA_parameter(double epsilon);
double weeks_chandler_ij(std::vector<double> position_i, std::vector<double> position_j);
double u_r(std::vector<double> r, std::vector<std::vector<double>> positions_ij);// the potential energy of a monomer
// at position r due to the excluded volume interaction (weeks chandler above) with existing monomers.


// select position for monomer from k trial positions
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions,
	std::vector<std::vector<double>>& existing_positions, double& weight, double& energy, int* selection_index=nullptr);

// overall probability of generating a set of monomers
double configuration_probability(std::vector<double>& energies, std::vector<double>& weights);

void rosenbluth_sample_helix_vectors(int n, std::vector<double>& origin, std::vector<double>& u, std::vector<double>& v, std::vector<std::vector<double>>& excluded_volume_positions, double& rosenbluth_factor, std::vector<std::vector<double>> existing_helix = {});


double rosenbluth_factor(std::vector<double>& weights, int k);

// very imporatnt growth functions. these are called many times to generate segments of teh chain according to a fixed
// end to end growth scheme or simply a random walk growth scheme. both use rosenbluth sampling.
std::vector<std::vector<double>> rosenbluth_grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l,
	std::vector<std::vector<double>>& excluded_volumes, std::vector<double>& energy_vals, std::vector<double>& weight_vals);

std::vector<std::vector<double>> rosenbluth_random_walk(std::vector<double>& starting_end, int monomers, std::vector<std::vector<double>>& excluded_volumes,
	std::vector<double>& energy_vals, std::vector<double>& weight_vals);

std::vector<std::vector<double>> grow_section(std::vector<std::vector<double>>& fixed_points, int segments, std::vector<std::vector<double>>& excluded_volumes, std::vector<double>& energy_vals, std::vector<double>& weight_vals, bool yamakawa, bool forward_move, std::vector<std::vector<double>> accepted_positions = {});

//*****************************************************************************************************************//
//*****************************************************************************************************************//
//*****************************************************************************************************************//
//*****************************************************************************************************************//
std::vector<double> rosenbluth_sample_translate(std::vector<std::vector<double>>& old_helix_positions, std::vector<std::vector<double>>& excluded_volume, double translation_distance, bool forward_move, double &rosenbluth_weight);

double rosenbluth_sample_corkscrew(helix_struct* dh, std::vector<std::vector<double>>& old_helix_positions, std::vector<std::vector<double>>& excluded_volume, double max_angle, bool forward_move, double &rosenbluth_weight);
std::vector<std::vector<double>> generate_corkscrew(helix_struct* dh, std::vector<std::vector<double>> helix_positions, double angle);

void rosenbluth_sample_spin(helix_struct* dh, std::vector<std::vector<double>>& old_helix_positions, std::vector<std::vector<double>>& excluded_volume, bool forward_move, double& rosenbluth_weight, double& selected_angle, std::vector<double>& selected_rot_axis);

std::vector<std::vector<double>> generate_spin(helix_struct* dh, std::vector<std::vector<double>>& helix_positions, double &angle, std::vector<double>& rotation_axis);

std::vector<double> centre_of_mass(std::vector<std::vector<double>>&helix_positions);

