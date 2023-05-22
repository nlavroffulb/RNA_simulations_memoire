#include <vector>
#include "math_functions.h"
#include "polymer_generation.h"
//const int trials{ 100 };


class rosenbluth_growth {
private:
	double P_R;// Rosenbluth probability distribution of generating the section
	double W_R;// Rosenbluth weight of the section
	double U;// the total energy of the section
	int l; // length of the section = number of monomers
	std::vector<std::vector<double>> positions;
	bool random_walk;//if random walk is true that means the section was grown 
	//with only one fixed point else it was grown between two fixed points
	int k{ 100 };
	std::vector<double> weights;
	std::vector<double> energies;
	std::vector<int> monomer_limits;



public:
	rosenbluth_growth() = default;
	rosenbluth_growth(double prob, double weight, std::vector<std::vector<double>> position_vectors, bool rw);
	rosenbluth_growth(double prob, double weight, double U, bool rw);
	rosenbluth_growth(std::vector<int>& limits, std::vector<double>& weight_vals, std::vector<double>& U_values, bool rw);
	~rosenbluth_growth();
	double get_config_prob();
	double get_rosen_weight();
	double get_energy();
	std::vector<std::vector<double>> get_positions();
	bool get_growth_method();

	double subsection_weight(std::vector<int> limit);
	double subsection_energy(std::vector<int> limit);
	double subsection_probability(std::vector<int> &limit);

	void modify_weights(std::vector<int> limits, std::vector<double> new_weights);
	void modify_energies(std::vector<int> limits, std::vector<double> new_energies);

};

double lennard_jones_ij(std::vector<double> position_i, std::vector<double> position_j);
double u_r(std::vector<double> r, std::vector<std::vector<double>> positions_ij);
double normalized_weight(std::vector<double> weights, int trials);// the normalized rosenbluth weight for a polymer

double selection_prob(std::vector<double>& trial_probabilities, double& trial_i);

// select position for monomer from k trial positions
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions,
	std::vector<std::vector<double>>& existing_positions, double& weight, double& energy);
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions,
	std::vector<std::vector<double>>& existing_positions);

// overall probability of generating a set of monomers
double configuration_probability(double total_energy, double normal_weight, int trials, int l);
double configuration_probability(std::vector<double>& energies, std::vector<double>& weights);
double configuration_probability(std::vector<std::vector<double>>& selected_positions, std::vector<double>& weights, std::vector<std::vector<double>> excluded_volume);

std::vector<std::vector<double>> select_helix(std::vector<std::vector<std::vector<double>>>& trial_helices,
	std::vector<std::vector<double>>& excluded_volume);


double rosenbluth_factor(std::vector<double> weights, int trials);


std::vector<std::vector<double>> rosenbluth_grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l,
	std::vector<std::vector<double>>& excluded_volumes, int trials, double& configuration_prob, rosenbluth_growth* section);
std::vector<std::vector<double>> rosenbluth_grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l,
	std::vector<std::vector<double>>& excluded_volumes, std::vector<double>& energy_vals, std::vector<double>& weight_vals);

std::vector<std::vector<double>> rosenbluth_random_walk(std::vector<double>& starting_end, int monomers, std::vector<std::vector<double>>& excluded_volumes, rosenbluth_growth* section);
std::vector<std::vector<double>> rosenbluth_random_walk(std::vector<double>& starting_end, int monomers, std::vector<std::vector<double>>& excluded_volumes,
	std::vector<double>& energy_vals, std::vector<double>& weight_vals);
