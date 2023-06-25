#include "rosenbluth_growth_class.h"
#include <algorithm>
//#include "rosenbluth_sampling.h"
double const a{ 1.0 };
int trials{ 5 };

rosenbluth_growth::rosenbluth_growth(int length)
{
	weights.resize(length);
	energies.resize(length);
}

rosenbluth_growth::rosenbluth_growth(std::vector<double>& weight_vals, std::vector<double>& energy_vals) {
	//monomer_limits = limits;
	weights = weight_vals;
	energies = energy_vals;
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
	if (weights.size() == 0) {
		return {};
	}
	else { return weights; }
}

void rosenbluth_growth::set_weights(std::vector<double> weight_vals)
{
	weights = weight_vals;
}




// returns the product of weights for a certain subset of monomers in the chain
double rosenbluth_growth::subsection_weight(std::vector<int> limit)
{
	if (limit[1] < limit[0]) {
		std::reverse(limit.begin(), limit.end());
	}

	std::vector<double> sub_weights = { weights.begin() + limit[0], weights.begin() + limit[1]+1 };
	//std::vector<double>   sub_weights(&weights[v_limits[0]], &weights[v_limits[1]]);
	return rosenbluth_factor(sub_weights, trials);
}

// returns the total energy for a certain subset of monomers in the chain.
double rosenbluth_growth::subsection_energy(std::vector<int> limit)
{
	if (limit[1] < limit[0]) {
		std::reverse(limit.begin(), limit.end());
	}

	std::vector<double> sub_energies = { weights.begin() + limit[0], weights.begin() + limit[1] };//might need to add 1 or -1 to these.
	//std::vector<double>   sub_weights(&weights[v_limits[0]], &weights[v_limits[1]]);
	return sum_of_elements(sub_energies);
}

// modify the weights member variable of the rosenbluth growth object for a certain limit.
void rosenbluth_growth::modify_weights(std::vector<int> limits, std::vector<double> new_weights)
{
	//for (auto i : new_weights) {
	//	if (isinf(i)) {
	//		std::cout << "stop" << std::endl;
	//	}
	//}
	int N = weights.size();
	int a{ limits[0] }, b{ limits[1] };

	if (limits[0] == limits[1] && new_weights.size()==1) {
		weights[limits[0]] = new_weights[0];
	}
	if (a > b) {
		std::reverse(new_weights.begin(), new_weights.end());
		a = limits[1];
		b = limits[0];
	}
	if (a < 0 || a > b) {
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
	if (limits[0] == limits[1] && new_energies.size() == 1) {
		energies[limits[0]] = new_energies[0];
	}

	if (a > b) {
		std::reverse(new_energies.begin(), new_energies.end());
		a = limits[1];
		b = limits[0];
	}

	if (a < 0 || a > b) {
		std::cout << "Invalid range specified" << std::endl;
		return;
	}
	for (int i = a; i <= b; i++) {
		energies[i] = new_energies[i - a]; // Modify the element in the range [a,b]
	}

}

rosenbluth_growth& rosenbluth_growth::operator=(const rosenbluth_growth& rhs)
{
	energies = rhs.energies;
	weights = rhs.weights;
	return *this;
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
		prod*=trials* exp(-energies[limit[0] + i]) / weights[limit[0] + i]; // the factor of k (trials) normalizes the weights.
		//prod *=  exp(-energies[limit[0] + i]) / weights[limit[0] + i]; // the factor of k (trials) normalizes the weights.

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

			}
		}
	}
	//std::cout << "energy " << summand << std::endl;
	return summand;
}


// this function is the trial section probability function. it takes the generating trial positions and the excluded volume as key 
// arguments. With those, it will select a trial and return the weight and energy of that trial (by reference, the weight and energy
// of the selected trial are arguments of the function).
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions,
	std::vector<std::vector<double>>& excluded_volume, double& weight, double& energy, int* selection_index)
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
		//std::cout << "trial " << j << std::endl;
		summand += exp(-u_r(trial_positions[j], excluded_volume));

		bins[j + 1] = summand;
	}
	weight = summand;
	int selection{ 0 };

	if (selection_index != nullptr) {// selection index is an optional argument. useful when calculating the k-1 trials of the old configuration
		selection = *selection_index;
	}
	else {
		double rn{ rand2(0,1) * summand };
		for (int k{ 1 }; k < bins.size(); k++) {
			if (rn < bins[k]) {
				selection = k - 1;
				break;
			}
		}
	}

	energy = u_r(trial_positions[selection], excluded_volume);
	return trial_positions[selection];
}

void rosenbluth_sample_helix_vectors(int n, std::vector<double> &origin, std::vector<double>& u, std::vector<double>& v, std::vector<std::vector<double>>& excluded_volume_positions, double &rosenbluth_factor, std::vector<std::vector<double>> existing_helix)
{
	std::vector<double> bins(trials + 1);

	bins[0] = 0;

	bool forward_move;
	existing_helix.size() == 0 ? forward_move = true : forward_move = false;

	int k{trials};
	std::vector<std::vector<double>> s_x(n), s_y(n), s(2*n);
	std::vector<double> boltzmann_weights(trials);
	std::vector<std::vector<double>> u_vectors(trials), v_vectors(trials);
	double rosenbluth_weight;
	double energy, cumulative_sum_of_weights{0};
	
	if (excluded_volume_positions.size() == 0) {
		sample_helix_vectors(u, v);
		rosenbluth_factor = 1;

	}
	else {
		if (forward_move == false) {
			k = k - 1;
			energy = 0;

			for (auto i : existing_helix) {
				energy = energy + u_r(i, excluded_volume_positions);

			}
			if (energy > 500) {
				boltzmann_weights.back() = 0.0;
			}
			else {
				boltzmann_weights.back() = exp(-energy);

			}

		}
		for (int i{ 0 }; i < k; i++) {
			energy = 0;
			sample_helix_vectors(u_vectors[i], v_vectors[i]);
			//dot_product(u_vectors[i], v_vectors[i]);
			generate(n, v_vectors[i], u_vectors[i], origin, s_x, s_y);
			s = s_x;
			s.insert(s.end(), s_y.begin(), s_y.end());
			//print_1d_doub_vec(s[1]);
			// calculate the energy and weights for each trial
			for (auto i : s) {
				energy = energy + u_r(i, excluded_volume_positions);
			}
			if (energy > 500) {
				boltzmann_weights[i] = 0.0;
			}
			else {
				boltzmann_weights[i] = exp(-energy);

			}
			cumulative_sum_of_weights = cumulative_sum_of_weights + boltzmann_weights[i];
			bins[i + 1] = cumulative_sum_of_weights;

		}
		int selection{ 0 };
		if (forward_move == true) {
			double rn{ rand2(0,1) * cumulative_sum_of_weights };
			for (int l{ 1 }; l < bins.size(); l++) {
				if (rn < bins[l]) {
					selection = l - 1;
					break;
				}
			}
			u = u_vectors[selection];
			v = v_vectors[selection];

		}
		else {
			selection = trials - 1;
		}
		//rosenbluth_factor = boltzmann_weights[selection]*trials/cumulative_sum_of_weights;
		rosenbluth_factor = cumulative_sum_of_weights / k;
	}
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
	std::vector<double> weights(l - 1); //the size of this vector is the number of monomers which are regrown (don't count fixed endpoints)
	std::vector<double> energies(l - 1); //the size of this vector is the number of monomers which are regrown (don't count fixed endpoints)


	//i is the monomer that is being grown. if you have 4 segments in total for example then you need
	//to add 3 monomers in addition to the endpoints so i=1,2,3. Note that including the endpoints we have
	//5 monomers in total in the flexible chain. 
	for (int i = 1; i <= l - 1; i++) {
		std::vector<std::vector<double>> trial_positions(trials);
		std::vector<double> trial_position(3);

		if (l - i > 1) {
			i_segments = l - i;
			for (int j{ 0 }; j < trials; j++) {
				trial_position = rejection_sample(old_position, ending_end, i_segments);
				trial_positions[j] = trial_position;
			}

			generated_position = select_trial(trial_positions, excluded_volumes, weights[i-1], energies[i-1]);
			excluded_volumes.push_back(generated_position);

			old_position = generated_position;
			positions[i] = generated_position;

		}
		else if (l - i == 1) {

			for (int k{ 0 }; k < trials; k++) {
				//std::cout << k << std::endl;
				trial_position = crankshaft_insertion(positions[l - 2], ending_end);
				trial_positions[k] = trial_position;
			}
			generated_position = select_trial(trial_positions, positions, weights[l-2], energies[l-2]);
			excluded_volumes.push_back(generated_position);

			positions[l - 1] = generated_position;

		}
	}
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
	std::vector<double> weights(monomers);
	std::vector<double> energies(monomers);


	std::vector<std::vector<double>> trial_positions(trials);
	for (int i{ 1 }; i <= monomers; i++) {
		for (int j{ 0 }; j < trials; j++) {
			sample_jump_direction(jump, 1);
			trial_positions[j] = vector_addition(old_position, jump);
		}
		generated_positions[i] = select_trial(trial_positions, excluded_volumes, weights[i-1], energies[i-1]);
		excluded_volumes.push_back(generated_positions[i]);
		old_position = generated_positions[i];

	}
	energy_vals = energies;
	weight_vals = weights;
	return generated_positions;
}

//the argument accepted positions must only be in the limit of growth we're considering. it should be a vector of size (segments + 1). actually need to be very careful about format of accepted positions vector, does it include the fixed points? is it in the right order ie if we're growing from alpha = 1.
std::vector<std::vector<double>> grow_section(std::vector<std::vector<double>>& fixed_points, int segments, std::vector<std::vector<double>>& excluded_volumes, std::vector<double>& energy_vals, std::vector<double>& weight_vals, bool yamakawa, bool forward_move, std::vector<std::vector<double>> accepted_positions)
	{
	std::vector<double> starting_end{ fixed_points[0] }, ending_end(3);

	std::vector<std::vector<double>> generated_positions(segments + 1);

	std::vector<double> generated_position(3);
	std::vector<double> old_position{ starting_end };

	std::vector<double> weights, energies;
	if (yamakawa == true) {
		generated_positions[0] = starting_end;
		weights.resize(segments - 1);//the size of this vector is the number of monomers which are regrown (don't count fixed endpoints)
		energies.resize(segments - 1);//the size of this vector is the number of monomers which are regrown (don't count fixed endpoints)
		ending_end = fixed_points[1];

	}
	else if (yamakawa == false) {
		generated_positions[0] = starting_end;
		weights.resize(segments);//the size of this vector is the number of monomers which are regrown (don't count fixed endpoints)
		energies.resize(segments);//the size of this vector is the number of monomers which are regrown (don't count fixed endpoints)

	}

	int k{ trials };
	// taking into account that we only generate k-1 trials for the old configuration.
	if (forward_move == false) { k = k - 1; }
	std::vector<std::vector<double>> trial_positions(trials);
	std::vector<double> trial_position(3);
	//std::cout << "checkpoint -1" << std::endl;

	int trial_index;
	// both ends fixed rosenbluth growth
	if (yamakawa == true) {
		int i_segments;
		//i is the monomer that is being grown. if you have 4 segments in total for example then you need
		//to add 3 monomers in addition to the endpoints so i=1,2,3. Note that including the endpoints we have
		//5 monomers in total in the flexible chain. 

		for (int i = 1; i <= segments - 1; i++) {

			if (segments - i > 1) {

				i_segments = segments - i;
				if (forward_move == false) {
					trial_positions[0] = accepted_positions[i - 1];

				}
				for (int j{ 0 }; j < k; j++) {
					trial_position = rejection_sample(old_position, ending_end, i_segments);
					forward_move == false ? trial_positions[j + 1] = trial_position : trial_positions[j] = trial_position;

				}

				if (forward_move == false) {
					int selection{ 0 };
					generated_position = accepted_positions[i - 1];
					select_trial(trial_positions, excluded_volumes, weights[i - 1], energies[i - 1], &selection);
				}
				else {
					generated_position = select_trial(trial_positions, excluded_volumes, weights[i - 1], energies[i - 1]);
				}
				excluded_volumes.push_back(generated_position);

				old_position = generated_position;
				generated_positions[i] = generated_position;

			}
			else if (segments - i == 1) {

				//if re generating the old configuration, we need to include the current positions in the trials.
				if (forward_move == false) {
					trial_positions[0] = accepted_positions[i - 1];
				}

				for (int j{ 0 }; j < k; j++) {
					trial_position = crankshaft_insertion(generated_positions[segments - 2], ending_end);
					forward_move == false ? trial_positions[j + 1] = trial_position : trial_positions[j] = trial_position;

				}
				if (forward_move == false) {
					int selection{ 0 };
					generated_position = accepted_positions[i - 1];

					//here we call select_trial without taking its return (a vector of the selected position). select trial also 
					//takes an extra, optional argument here which is the selection. We force it to select the current position as
					// the trial. all we're interested in when we call select trial here is that it returns the weight and energy
					// by reference.
					select_trial(trial_positions, excluded_volumes, weights[i - 1], energies[i - 1], &selection);
				}
				else {
					generated_position = select_trial(trial_positions, excluded_volumes, weights[i - 1], energies[i - 1]);
				}

				excluded_volumes.push_back(generated_position);

				generated_positions[segments - 1] = generated_position;

			}
		}
		generated_positions[segments] = ending_end;



	}
	// random walk rosenbluth growth
	else if (yamakawa == false) {
		std::vector<double> jump(3);


		for (int i{ 1 }; i <= segments; i++) {

			if (forward_move == false) {
				trial_positions[0] = accepted_positions[i - 1];
			}

			for (int j{ 0 }; j < k; j++) {
				sample_jump_direction(jump, 1);
				forward_move == false ? trial_positions[j + 1] = vector_addition(old_position, jump) : trial_positions[j] = vector_addition(old_position, jump);
			}
			if (forward_move == false) {
				int selection{ 0 };
				generated_position = accepted_positions[i - 1];
				//here we call select_trial without taking its return (a vector of the selected position). select trial also 
				//takes an extra, optional argument here which is the selection. We force it to select the current position as
				// the trial. all we're interested in when we call select trial here is that it returns the weight and energy
				// by reference.
				select_trial(trial_positions, excluded_volumes, weights[i - 1], energies[i - 1], &selection);
			}
			else {
				generated_position = select_trial(trial_positions, excluded_volumes, weights[i - 1], energies[i - 1]);
			}
			generated_positions[i] = generated_position;

			excluded_volumes.push_back(generated_positions[i]);
			old_position = generated_positions[i];

		}

	}

	energy_vals = energies;
	weight_vals = weights;
	return generated_positions;
}


std::vector<double> rosenbluth_sample_translate(std::vector<std::vector<double>>& old_helix_positions, std::vector<std::vector<double>>& excluded_volume, double translation_distance, bool forward_move, double &rosenbluth_weight)
{
	int k{trials};
	std::vector<double> trial_translation(3);
	std::vector<std::vector<double>> trial_helix_positions(old_helix_positions.size());
	double cumulative_weight{ 0 }, trial_energy{ 0 };
	if (forward_move == false) {
		k = trials-1;
		for (auto i : old_helix_positions) {
			trial_energy += u_r(i, excluded_volume);
			if (trial_energy > 100) {//the weight will be 0 in this case. so we will catch this later in the acceptance rule and reject the move.
				rosenbluth_weight = 0;
				return{ 0,0,0 };
			}
			cumulative_weight = exp(-trial_energy);
		}
	}
	std::vector<std::vector<double>> translations(k);

	std::vector<double> bins(k + 1);
	bins[0] = 0;

	for (int i = 0; i < k; i++) {
		if (i == k - 2) {
			std::cout << "stop" << std::endl;
		}

		trial_energy = 0;
		sample_jump_direction(trial_translation, translation_distance);
		translations[i] = trial_translation;
		for (int j{ 0 }; j < old_helix_positions.size();j++){
			trial_helix_positions[i] = vector_addition(old_helix_positions[i], trial_translation);
			trial_energy += u_r(trial_helix_positions[i], excluded_volume);
		}
		if (trial_energy > 200) {//the weight will be 0 in this case. so we will catch this later in the acceptance rule and reject the move.
			rosenbluth_weight = 0;
			return {0,0,0};
		}else{
			cumulative_weight += exp(-trial_energy);
			bins[i + 1] = cumulative_weight;
		}
		
	}
	int selection;
	if (forward_move == true) {
		double rn{ rand2(0,1) * cumulative_weight };
		for (int l{ 1 }; l < bins.size(); l++) {
			if (rn < bins[l]) {
				selection = l - 1;
				break;
			}
		}
		rosenbluth_weight = cumulative_weight;
		return translations[selection];
		
	}
	else {
		rosenbluth_weight = cumulative_weight;
		return { 0,0,0 };
	}
}

double rosenbluth_sample_corkscrew(helix_struct* dh, std::vector<std::vector<double>>& old_helix_positions, std::vector<std::vector<double>>& excluded_volume, double max_angle, bool forward_move, double& rosenbluth_weight)
{
	int k{ trials };
	double trial_angle;
	std::vector < std::vector<double>> trial_helix_positions(old_helix_positions.size());
	double cumulative_weight{ 0 }, trial_energy{ 0 };
	if (forward_move == false) {
		k = trials - 1;
		for (auto i : old_helix_positions) {
			trial_energy += u_r(i, excluded_volume);
			if (trial_energy > 200) {//the weight will be 0 in this case. so we will catch this later in the acceptance rule and reject the move.
				rosenbluth_weight = 0;
				return 0.0;
			}
			cumulative_weight = exp(-trial_energy);
		}
	}
	std::vector<double> trial_angles(k);
	std::vector<double> bins(k + 1);
	bins[0] = 0;


	for (int i = 0; i < k; i++)
	{
		if (i == k - 2) {
			std::cout << "stop" << std::endl;
		}
		trial_energy = 0;
		trial_angle = rand2(0, max_angle);
		trial_helix_positions = generate_corkscrew(dh, old_helix_positions, trial_angle);
		trial_angles[i] = trial_angle;

		for (int j{ 0 }; j < trial_helix_positions.size(); j++) {
			trial_energy += u_r(trial_helix_positions[i], excluded_volume);
		}
		if (trial_energy > 200) {//the weight will be 0 in this case. so we will catch this later in the acceptance rule and reject the move.
			rosenbluth_weight = 0;
			return 0.0;
		}
		else {
			cumulative_weight += exp(-trial_energy);
			bins[i + 1] = cumulative_weight;
		}

	}
	int selection;
	if (forward_move == true) {
		double rn{ rand2(0,1) * cumulative_weight };
		for (int l{ 1 }; l < bins.size(); l++) {
			if (rn < bins[l]) {
				selection = l - 1;
				break;
			}
		}
		rosenbluth_weight = cumulative_weight;
		return trial_angles[selection];

	}
	else {
		rosenbluth_weight = cumulative_weight;
		return 0.0;
	}

}

std::vector<std::vector<double>> generate_corkscrew(helix_struct* dh, std::vector<std::vector<double>> helix_positions, double angle)
{
	std::vector<std::vector<double>> rotated_positions(helix_positions.size());

	std::vector<int> s{ dh->get_monomers() };
	std::vector<double> v{ dh->get_v() };
	normalize(v);

	int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };
	std::vector<double> running_centre(3);
	beta == 0 ? running_centre = dh->get_rc(beta) : running_centre = dh->get_rc(1 - beta);

	int n{ abs(s[0] - s[1]) + 1 };

	// starting pair to be rotated depends on beta
	int m1, m2;
	if (beta == 0)
	{
		m1 = 0;
		m2 = helix_positions.size()-1;
	}
	else {
		m1 = helix_positions.size()/2;
		m2 = -1+ helix_positions.size() / 2;
	}
	int gamma{ static_cast<int>(pow(-1,beta)) };
	double pitch{ 0.28 };
	std::vector<double> temp{ multiplication_by_scalar(pitch, v) };

	std::vector<double> p1{ helix_positions[m1] }, p2{ helix_positions[m2] };

	p1 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);
	p2 = rotate_by_angle_about_general_axis(p2, v, angle, running_centre);
	rotated_positions[m1] = p1;
	rotated_positions[m2] = p2;
	//chain[m1]->change_position(p1), chain[m2]->change_position(p2);
	for (int i = 1; i < n; i++)
	{
		running_centre = vector_addition(running_centre, temp);

		m1 = m1 + gamma;
		m2 = m2 - gamma;

		p1 = helix_positions[m1];
		p2 = helix_positions[m2];

		p1 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);
		p2 = rotate_by_angle_about_general_axis(p2, v, angle, running_centre);

		rotated_positions[m1] = p1, rotated_positions[m2] = p2;

	}
	return rotated_positions;
}

void rosenbluth_sample_spin(helix_struct* dh, std::vector<std::vector<double>>& old_helix_positions, std::vector<std::vector<double>>& excluded_volume, bool forward_move, double& rosenbluth_weight, double& selected_angle, std::vector<double>& selected_rot_axis)
{
	int k{ trials };
	std::vector < std::vector<double>> trial_helix_positions(old_helix_positions.size());
	double cumulative_weight{ 0 }, trial_energy{ 0 };
	if (forward_move == false) {
		k = trials - 1;
		for (auto i : old_helix_positions) {
			trial_energy += u_r(i, excluded_volume);
			if (trial_energy > 200) {//the weight will be 0 in this case. so we will catch this later in the acceptance rule and reject the move.
				rosenbluth_weight = 0;
				return;
			}
			cumulative_weight = exp(-trial_energy);
		}
	}
	double trial_angle;
	std::vector<double> trial_rotation_axis(3);
	std::vector<double> trial_angles(k);
	std::vector<std::vector<double>> trial_rotation_axes(k);
	std::vector<double> bins(k + 1);
	bins[0] = 0;

	for (int i = 0; i < k; i++) {
		if (i == k - 2) {
			std::cout << "stop" << std::endl;
		}

		trial_energy = 0;
		trial_helix_positions = generate_spin(dh, old_helix_positions, trial_angle, trial_rotation_axis);

		trial_angles[i] = trial_angle;
		trial_rotation_axes[i] = trial_rotation_axis;

		for (int j{ 0 }; j < old_helix_positions.size(); j++) {
			trial_energy += u_r(trial_helix_positions[i], excluded_volume);
		}
		if (trial_energy > 200) {//the weight will be 0 in this case. so we will catch this later in the acceptance rule and reject the move.
			rosenbluth_weight = 0;
			return;
		}
		else {
			cumulative_weight += exp(-trial_energy);
			bins[i + 1] = cumulative_weight;
		}


	}

	int selection;
	if (forward_move == true) {
		double rn{ rand2(0,1) * cumulative_weight };
		for (int l{ 1 }; l < bins.size(); l++) {
			if (rn < bins[l]) {
				selection = l - 1;
				break;
			}
		}
		rosenbluth_weight = cumulative_weight;
		selected_angle = trial_angles[selection];
		selected_rot_axis = trial_rotation_axes[selection];
		return;

	}
	else {
		rosenbluth_weight = cumulative_weight;
		return;
	}

}

std::vector<std::vector<double>> generate_spin(helix_struct* dh, std::vector<std::vector<double>>& helix_positions, double &angle, std::vector<double>& rotation_axis)
{

	std::vector<int> s{ dh->get_monomers() };

	std::vector<double> com{ centre_of_mass(helix_positions)};

	std::vector<double> v_p{ dh->get_v() };
	double theta{ rand2(0,4 * atan(1)) }; //max rotation is pi.
	std::vector<double> rot_axis{ sample_unit_perp_vector(v_p) };
	rotation_axis = rot_axis;

	std::vector<std::vector<double>> rotated_positions(helix_positions.size());
	std::vector<double> temp(3);
	for (int i{ 0 }; i < helix_positions.size(); i++) {
		temp = helix_positions[i];
		rotated_positions[i] = rotate_by_angle_about_general_axis(temp, rot_axis, theta, com);
	}

	rotation_axis = rot_axis;
	angle = theta;

	return rotated_positions;
	
}

std::vector<double> centre_of_mass(std::vector<std::vector<double>>& helix_positions)
{
	std::vector<double> r_com{ helix_positions[0]};

	std::vector<double> temp(3);
	for (int i{ 1 }; i < helix_positions.size(); i++) {
		temp = helix_positions[i];
		r_com = vector_addition(r_com, temp);
	}
	int k{ static_cast<int>(helix_positions.size()) };
	std::for_each(r_com.begin(), r_com.end(), [k](double& c) { c /= k; });

	return r_com;
}


double rosenbluth_factor(std::vector<double>& weights, int trials)
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
