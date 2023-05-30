#include "rosenbluth_growth_class.h"
//#include "rosenbluth_sampling.h"
double const a{ 1.0 };
int trials{ 100 };

rosenbluth_growth::rosenbluth_growth(std::vector<double>& weight_vals, std::vector<double>& energy_vals) {
	//monomer_limits = limits;
	weights = weight_vals;
	energies = energy_vals;
	//random_walk = rw;
	l = weight_vals.size();
}
rosenbluth_growth::~rosenbluth_growth() {
	std::cout << "Unbound section destructor called" << std::endl;
}

double rosenbluth_growth::get_config_prob()
{
	return configuration_probability(energies, weights);
}

double rosenbluth_growth::get_rosen_weight()
{
	return rosenbluth_factor(weights, trials);
}

double rosenbluth_growth::get_energy()
{
	return sum_of_elements(energies);
}

double rosenbluth_growth::get_Ui(int index)
{
	return energies[index];
}

double rosenbluth_growth::get_Wi(int index)
{
	return weights[index];
}

std::vector<double> rosenbluth_growth::get_energies()
{
	return energies;
}

void rosenbluth_growth::set_energies(std::vector<double> energy_vals)
{
	energies = energy_vals;
	return;
}

std::vector<double> rosenbluth_growth::get_weights()
{
	return weights;
}

void rosenbluth_growth::set_weights(std::vector<double> weight_vals)
{
	weights = weight_vals;
}



bool rosenbluth_growth::get_growth_method()
{
	return random_walk;
}

// returns the product of weights for a certain subset of monomers in the chain
double rosenbluth_growth::subsection_weight(std::vector<int> limit)
{

	std::vector<double> sub_weights = { weights.begin() + limit[0], weights.begin() + limit[1] };
	//std::vector<double>   sub_weights(&weights[v_limits[0]], &weights[v_limits[1]]);
	return rosenbluth_factor(sub_weights, trials);
}

// returns the total energy for a certain subset of monomers in the chain.
double rosenbluth_growth::subsection_energy(std::vector<int> limit)
{
	std::vector<double> sub_energies = { weights.begin() + limit[0], weights.begin() + limit[1] };//might need to add 1 or -1 to these.
	//std::vector<double>   sub_weights(&weights[v_limits[0]], &weights[v_limits[1]]);
	return sum_of_elements(sub_energies);
}

// modify the weights member variable of the rosenbluth growth object for a certain limit.
void rosenbluth_growth::modify_weights(std::vector<int> limits, std::vector<double> new_weights)
{
	int N = weights.size();
	int a{ limits[0] }, b{ limits[1] };
	if (a > b) {
		std::reverse(new_weights.begin(), new_weights.end());
		a = limits[1];
		b = limits[0];
	}
	if (a < 0 || b >= N || a > b) {
		std::cout << "Invalid range specified" << std::endl;
		return;
	}
	for (int i = a; i <= b; i++) {
		weights[i] =new_weights[i-a]; // Modify the element in the range [a,b]
	}
}


void rosenbluth_growth::modify_energies(std::vector<int> limits, std::vector<double> new_energies)
{
	int N = energies.size();
	int a{ limits[0] }, b{ limits[1] };
	if (a > b) {
		std::reverse(new_energies.begin(), new_energies.end());
		a = limits[1];
		b = limits[0];
	}

	if (a < 0 || b >= N || a > b) {
		std::cout << "Invalid range specified" << std::endl;
		return;
	}
	for (int i = a; i <= b; i++) {
		energies[i] = new_energies[i - a]; // Modify the element in the range [a,b]
	}

}

// calculate the probability of generating a certain section of the chain, worth noting that this is dependent on the order
// that the chain was grown although this is not explicitly taken into account in this function. (it should have been configured 
// beforehand in such a way that it is taken into account).
double rosenbluth_growth::subsection_probability(std::vector<int> &limit) {
	//int l{ abs(limit[1] - limit[0]) + 1 };
	int l{ abs(limit[1] - limit[0])+1};

	double prod{ 1 };

	if (limit[1] < limit[0]) {
		std::reverse(limit.begin(), limit.end());
	}
	//assuming the limit is from high to low.
	for (int i{ 0 }; i < l; i++) {
		prod*=exp(-energies[limit[0] + i]) / weights[limit[0] + i];
	}
	if (limit[1] < limit[0]) {
		std::reverse(limit.begin(), limit.end());
	}
	return prod;
}

// this is a function which calculates the rosenbluth distribution probability of a chain given the weights
// and energies of that chain.
double configuration_probability(std::vector<double>& energies, std::vector<double>& weights) {
	double product{ 1 };

	for (int i{ 0 }; i < energies.size(); i++) {
		product *= exp(-energies[i]) / weights[i];
	}
	return product;
}


//*****************************************************************************************************************//
//*****************************************************************************************************************//
//*****************************************************************************************************************//
//*****************************************************************************************************************//

double weeks_chandler_ij(std::vector<double> position_i, std::vector<double> position_j) {// this is actually the 
	// weeks chandler andersen intermolecular pair potential (which is a function of the LJ potential)
	if (!position_i.empty() && !position_j.empty()) {
		double r_ij{ dist_2_points3d(position_i,position_j) };

		if (r_ij <= a * pow(2, 1.0 / 6.0)) {
			//if (4 * (pow(a / r_ij, 12) - pow(a / r_ij, 6)) + 1 > 1000) {
			//	std::cout << "stop" << std::endl;
			//}

			return 4 * (pow(a / r_ij, 12) - pow(a / r_ij, 6)) + 1;
		}
		else if (r_ij > a * pow(2, 1.0 / 6.0)) {
			return 0;
		}

	}
}

double u_r(std::vector<double> r, std::vector<std::vector<double>> positions_ij) {
	double summand{ 0 };
	for (int j{ 0 }; j < positions_ij.size(); j++) {
		//sometimes we will calculate u(r) and positions_ij will be basically empty.
		//issue is that we declare its size at the start, so it can still loop through the elements
		//even if the elements themselves are empty.

		if (positions_ij[j] == r) {// we don't count ev with itself.
			continue;
		}
		else {
			//if (dist_2_points3d(positions_ij[j], r)) {

			//}
			if (!positions_ij[j].empty() && !r.empty()) {
				summand += weeks_chandler_ij(r, positions_ij[j]);
				//if (summand > 1000) {
				//	std::cout << "stop" << std::endl;
				//}
			}
		}
	}
	return summand;
}


// this function is the trial section probability function. it takes the generating trial positions and the excluded volume as key 
// arguments. With those, it will select a trial and return the weight and energy of that trial (by reference, the weight and energy
// of the selected trial are arguments of the function).
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions,
	std::vector<std::vector<double>>& excluded_volume, double& weight, double& energy)
{
	std::vector<double> bins(trial_positions.size() + 1);

	bins[0] = 0;
	double summand{ 0 };
	for (auto i : excluded_volume) {
		if (i.size() == 0) {
			break;
		}
	}
	for (int j{ 0 }; j < trial_positions.size(); j++) {
		summand += exp(-u_r(trial_positions[j], excluded_volume));

		bins[j + 1] = summand;
	}
	if (summand == 0) {
		std::cout << "stop " << std::endl;
		double summand{ 0 };

		for (int j{ 0 }; j < trial_positions.size(); j++) {
			summand += exp(-u_r(trial_positions[j], excluded_volume));

			bins[j + 1] = summand;
		}

	}

	weight = summand;
	double rn{ rand2(0,1) * summand };
	int selection{ 0 };
	for (int k{ 1 }; k < bins.size(); k++) {
		if (rn < bins[k]) {
			selection = k - 1;
			break;
		}
	}
	energy = u_r(trial_positions[selection], excluded_volume);
	return trial_positions[selection];
}

// this is a function that hasn't been revised in a while but in principal it would apply a Rosenbluth sampling
//scheme to the growth of a helix. we grow several helixes and calculate the excluded volume interaction energy of 
//each of them then do a typical roulette method to select.
std::vector<std::vector<double>> select_helix(std::vector<std::vector<std::vector<double>>>& helices, std::vector<std::vector<double>>& existing_positions)
{

	std::vector<double> bins(helices.size() + 1);

	bins[0] = 0;
	double summand{ 0 };
	for (int i{ 0 }; i < helices.size(); i++) {
		for (int j{ 0 }; j < helices[i].size(); j++) {
			summand += exp(-u_r(helices[i][j], existing_positions));
		}
		bins[i + 1] = summand;
	}
	double rn{ rand2(0,1) * summand };
	int selection{ 0 };
	for (int k{ 1 }; k < bins.size(); k++) {
		if (rn < bins[k]) {
			selection = k - 1;
			break;
		}
	}
	return helices[selection];
}

// the function which applies a Rosenbluth sampled, fixed end to end growth of monomers. the function returns the positions of the 
// monomers and their Rosenbluth weights and energies. The positions are returned in the standard way while the energies and weights 
// are returned by reference, they are arguments of the function.
std::vector<std::vector<double>> rosenbluth_grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l,
	std::vector<std::vector<double>>& excluded_volumes, std::vector<double>& energy_vals, std::vector<double>& weight_vals) {

	std::vector<std::vector<double>> positions(l + 1);
	std::vector<double> generated_position(3);
	std::vector<double> old_position{ starting_end };
	positions[0] = starting_end;
	int i_segments;
	std::vector<double> weights(l + 1);
	weights[0] = 1;
	double normal_weight;
	std::vector<double> energies(l + 1);
	energies[0] = 0;

	//i is the monomer that is being grown. if you have 4 segments in total for example then you need
	//to add 3 monomers in addition to the endpoints so i=1,2,3. Note that including the endpoints we have
	//5 monomers in total in the flexible chain. 
	for (int i = 1; i <= l - 1; i++) {
		std::vector<std::vector<double>> trial_positions(trials);
		std::vector<double> trial_position(3);

		if (l - i > 1) {
			i_segments = l - i;
			//if(i!=1)
			if (i == 1) {
				trial_position = rejection_sample(old_position, ending_end, i_segments);
				energies[1] = u_r(trial_position, excluded_volumes);
				weights[1] = trials * exp(-energies[0]);
				excluded_volumes.push_back(trial_position);
				generated_position = trial_position;
			}
			else {
				for (int j{ 0 }; j < trials; j++) {
					trial_position = rejection_sample(old_position, ending_end, i_segments);
					trial_positions[j] = trial_position;
				}
				
				generated_position = select_trial(trial_positions, excluded_volumes, weights[i], energies[i]);
				if (weights[i] == 0) {
					std::cout << "problem" << std::endl;
				}
				excluded_volumes.push_back(generated_position);
			}
			old_position = generated_position;
			positions[i] = generated_position;

		}
		else if (l - i == 1) {

			for (int k{ 0 }; k < trials; k++) {
				//std::cout << k << std::endl;
				trial_position = crankshaft_insertion(positions[l - 2], ending_end);
				trial_positions[k] = trial_position;
			}
			generated_position = select_trial(trial_positions, positions, weights[l-1], energies[l-1]);
			excluded_volumes.push_back(generated_position);

			positions[l - 1] = generated_position;

		}
	}
	weights.back() = 1;
	energies.back() = 0;
	positions[l] = ending_end;
	energy_vals = energies;
	weight_vals = weights;

	return positions;
}


// the function which applies a Rosenbluth sampled, random walk growth of monomers. the function returns the positions of the 
// monomers and their Rosenbluth weights and energies. The positions are returned in the standard way while the energies and weights 
// are returned by reference, they are arguments of the function.
std::vector<std::vector<double>> rosenbluth_random_walk(std::vector<double>& starting_end, int monomers, std::vector<std::vector<double>>& excluded_volumes,
	std::vector<double>& energy_vals, std::vector<double>& weight_vals)
{
	std::vector<std::vector<double>> generated_positions(monomers + 1);//plus one because we're going to include
	// the starting point in the return vector of positions
	generated_positions[0] = starting_end;
	std::vector<double> old_position{ starting_end };
	std::vector<double> jump(3);
	std::vector<double> weights(monomers + 1);
	std::vector<double> energies(monomers + 1);

	energies[0] = 1;
	weights[0] = 1;
	excluded_volumes.push_back(starting_end);


	std::vector<std::vector<double>> trial_positions(trials);
	for (int i{ 1 }; i <= monomers; i++) {
		for (int j{ 0 }; j < trials; j++) {
			sample_jump_direction(jump, 1);
			trial_positions[j] = vector_addition(old_position, jump);
		}
		generated_positions[i] = select_trial(trial_positions, excluded_volumes, weights[i], energies[i]);
		excluded_volumes.push_back(generated_positions[i]);
		old_position = generated_positions[i];

	}
	energy_vals = energies;
	weight_vals = weights;
	return generated_positions;
}


double rosenbluth_factor(std::vector<double> weights, int trials)
{
	double prod{ 1 };
	int l{ static_cast<int>(weights.size()) };

	for (int i{ 0 }; i < l; i++) {
		prod *= weights[i] / trials;
	}
	return prod;
}


//*****************************************************************************************************************//
//*****************************************************************************************************************//
//*****************************************************************************************************************//
//*****************************************************************************************************************//
