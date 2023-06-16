#include <vector>
#include "math_functions.h"
#include "polymer_generation.h"
#include "helix_generation.h"

class rosenbluth_growth {
private:
	int l; // length of the section = number of monomers
	bool random_walk;//if random walk is true that means the section was grown 
	//with only one fixed point else it was grown between two fixed points
	std::vector<double> weights;
	std::vector<double> energies;
	//std::vector<int> monomer_limits;



public:

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
	bool get_growth_method();

	double subsection_weight(std::vector<int> limit);
	double subsection_energy(std::vector<int> limit);
	double subsection_probability(std::vector<int> &limit);

	void modify_weights(std::vector<int> limits, std::vector<double> new_weights);
	void modify_energies(std::vector<int> limits, std::vector<double> new_energies);

	rosenbluth_growth& operator=(const rosenbluth_growth& rhs);
};

double weeks_chandler_ij(std::vector<double> position_i, std::vector<double> position_j);
double u_r(std::vector<double> r, std::vector<std::vector<double>> positions_ij);// the potential energy of a monomer
// at position r due to the excluded volume interaction (weeks chandler above) with existing monomers.


// select position for monomer from k trial positions
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions,
	std::vector<std::vector<double>>& existing_positions, double& weight, double& energy);

// overall probability of generating a set of monomers
double configuration_probability(std::vector<double>& energies, std::vector<double>& weights);

void rosenbluth_sample_helix(int n, std::vector<double>& origin, std::vector<double>& u, std::vector<double>& v, std::vector<std::vector<double>> excluded_volume_positions);

// hasn't been implemented.
std::vector<std::vector<double>> select_helix(std::vector<std::vector<std::vector<double>>>& trial_helices,
	std::vector<std::vector<double>>& excluded_volume);


double rosenbluth_factor(std::vector<double> weights, int k);

// very imporatnt growth functions. these are called many times to generate segments of teh chain according to a fixed
// end to end growth scheme or simply a random walk growth scheme. both use rosenbluth sampling.
std::vector<std::vector<double>> rosenbluth_grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l,
	std::vector<std::vector<double>>& excluded_volumes, std::vector<double>& energy_vals, std::vector<double>& weight_vals);

std::vector<std::vector<double>> rosenbluth_random_walk(std::vector<double>& starting_end, int monomers, std::vector<std::vector<double>>& excluded_volumes,
	std::vector<double>& energy_vals, std::vector<double>& weight_vals);
