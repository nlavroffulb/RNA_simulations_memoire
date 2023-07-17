#include "polymer_generation.h"
//#include "rosenbluth_sampling.h"

// l = number of segments to be regrown
// i ranges from 1 to l
// the total number of monomers is l+1 (counting the start and end monomers)

double const a{ 1.0 };



double ideal_chain_pdf(double &e2e_distance, int n_segments)//'yamakawa' function.
// protected against factorial overflow if n>175.
{
	if (n_segments < 170) {
		double sum_value{ 0 };
		for (int k{ 0 }; k <= (n_segments - (e2e_distance / a)) / 2; k++) {
			sum_value += std::pow(-1, k) * choose(n_segments, k) * std::pow(n_segments - 2.0 * k - (e2e_distance / a), n_segments - 2);
		}
		return sum_value / (std::pow(2, n_segments + 1) * factorial(n_segments - 2) * pi * std::pow(a, 2) * e2e_distance);

	}
	else if (n_segments > 170) {
		return std::pow(1.5 / (pi * n_segments), 3.0 / 2.0) * exp(-1.5 * std::pow(e2e_distance, 2) / n_segments);
	}
}

double segment_grow_prob(double new_e2e, double old_e2e, int l_minus_i) {
	return ideal_chain_pdf(new_e2e, l_minus_i) / (4*pi*ideal_chain_pdf(old_e2e, l_minus_i +1));
	//return ideal_chain_pdf(new_e2e, l_minus_i);

}
double unnorm_segment_grow_prob(double new_e2e, int l_minus_i) {
	return ideal_chain_pdf(new_e2e, l_minus_i);

}
bool overstretch(int l_minus_i, double &trial_e2e) {
	if (trial_e2e > (l_minus_i)* a) {
		return true;
	}
	else { return false; }
}

double sample_e2e(double init_e2e_distance) {
	return init_e2e_distance + rand2(-a, a);
}

void sample_jump_direction(std::vector<double> &jump, double length) {
	length = a;

	double phi{ rand2(0,2 * pi) }, theta{ acos(1-2*rand2(0,1)) };//correct sampling of cos theta using 
	//inverse transform sampling

	jump = { a * sin(theta) * cos(phi),a * sin(theta) * sin(phi),a * cos(theta) };

}

std::vector<double> rejection_sample(std::vector<double> &initial, std::vector<double> &N_position, int segments_to_regrow) {

	if (dist_2_points3d(N_position, initial) > segments_to_regrow+1) {
		std::cout << "stop rejec sample" << std::endl;
	}
	if (segments_to_regrow > 1) {
		//pdf max
		std::vector<double> temp{ vector_subtraction(N_position,initial) };
		double init_e2e{ vector_modulus(temp) };
		double pdf_max{ unnorm_segment_grow_prob(std::abs(init_e2e - a),segments_to_regrow) };

		//trial jump
		std::vector<double> jump(3),test(3);
		sample_jump_direction(jump, a);
		std::vector<double> trial_position{ vector_addition(initial,jump) };
		temp = vector_subtraction(N_position, trial_position);
		double trial_e2e{ vector_modulus(temp) };


		int over_count{ 0 };
		double Y{ rand2(0,pdf_max) };

		//while (overstretch(segments_to_regrow-1, trial_e2e) == true || Y > segment_grow_prob(trial_e2e, init_e2e, segments_to_regrow)|| (std::abs(segments_to_regrow-1-trial_e2e)<0.1 && segments_to_regrow>3)) {
			//std::cout << "rejected" << std::endl;
		while (overstretch(segments_to_regrow, trial_e2e) == true || Y > unnorm_segment_grow_prob(trial_e2e, segments_to_regrow)) {

			Y = rand2(0, pdf_max);

			sample_jump_direction(jump, a);
			trial_position = vector_addition(initial, jump);
			test = vector_subtraction(N_position, initial);
			double init_e2e{ vector_modulus(test) };
			temp = vector_subtraction(N_position, trial_position);
			trial_e2e = vector_modulus(temp);
			over_count++;
			if (over_count == 10000) {
				std::cout << "might be stuck" << std::endl;
			}
			if (over_count == 100000000) {
				double_vector_to_txt("Results/stuck_in_rejection_sample.txt", {});
				over_count++;
				std::cout << "stuck" << std::endl;
				if (segments_to_regrow - 1 > trial_e2e) {
					std::cout << "outie" << std::endl;
				}
			}
		}//std::cout << "rejections" << over_count << std::endl;
		//if (over_count > 10000) {
		//	over_count++;
		//	std::cout << "ESCAPE AFTER STUCK" << std::endl;
		//}
		if (dist_2_points3d(N_position, trial_position) > segments_to_regrow) {
			std::cout << "stop rejec sample" << std::endl;
		}

		return trial_position;

	}
	else {
		std::cout << "The last segment needs to be regrown by crankshaft insertion" << std::endl;
		
	}

}

//this could break if we have an r_e2e vector with 0 components.
std::vector<double> crankshaft_insertion(std::vector<double> &start_pos, std::vector<double> &N_pos) {
	std::vector<double> r_e2e{ vector_subtraction(N_pos,start_pos) };

	double mod_r_e2e{ vector_modulus(r_e2e) };
	double r_c{ sqrt(a * a - std::pow(mod_r_e2e / 2,2)) };
	r_e2e = normalize(r_e2e);
	double theta = rand2(0, 2 * pi);
	double phi = acos(mod_r_e2e / (2 * a));
	double x{ cos(theta) }, y{ sin(theta) };

	double z{ -1 / r_e2e[2] * (x * r_e2e[0] + y * r_e2e[1]) };

	std::vector<double> v_c{ x,y,z };
	v_c = normalize(v_c);

	r_e2e = multiplication_by_scalar(mod_r_e2e / 2, r_e2e);
	v_c = multiplication_by_scalar(r_c, v_c);

	std::vector<double> r_Nminus1{ vector_addition(r_e2e,v_c) };

	r_Nminus1 = vector_addition(start_pos, r_Nminus1);
	//for (auto i : r_Nminus1) {
	//	if (isnan(i)) {
	//		std::cout << "stop" << std::endl;
	//	}
	//}
	return r_Nminus1;
};


std::vector<std::vector<double>> grow_chain(std::vector<double>& starting_end, 
	std::vector<double>& ending_end, int l, 
	std::vector<std::vector<double>>& excluded_volumes) {

	std::vector<std::vector<double>> positions(l+1);
	std::vector<double> generated_position(3);
	std::vector<double> old_position{ starting_end };
	positions[0]=starting_end;
	int i_segments;

	
	//i is the monomer that is being grown. if you have 4 segments in total for example then you need
	//to add 3 monomers in addition to the endpoints so i=1,2,3. Note that including the endpoints we have
	//5 monomers in total in the flexible chain. 
	for (int i=1; i <= l-1; i++) {
		std::vector<double> trial_position(3);

		//std::cout <<"regrowth stage " << i << std::endl;
		if (l - i > 1) {
			i_segments = l - i;
			//if(i!=1)
			if (excluded_volumes.size() > 0) {
				trial_position = rejection_sample(old_position, ending_end, i_segments);

				//excluded_volumes.push_back()
				//generated_position = select_trial(trial_positions, excluded_volumes);
				generated_position = trial_position;

				excluded_volumes.push_back(generated_position);

			}
			else { generated_position = rejection_sample(old_position, ending_end, i_segments); }
			old_position = generated_position;
			positions[i] = generated_position;

		}
		else if(l-i==1) {
			//std::cout << "Crankshaft regime" << std::endl;
			trial_position = crankshaft_insertion(positions[l - 2], ending_end);
			generated_position = trial_position;

			excluded_volumes.push_back(generated_position);

			positions[l-1] = generated_position;

		}
	}
	positions[l] = ending_end;

	//for (int i = 0; i < positions.size(); i++)

	//{
	//	std::cout << i;

	//	std::cout << "{";
	//	for (int j = 0; j < positions[i].size(); j++)
	//	{
	//		std::cout << positions[i][j] << " ";
	//	}
	//	std::cout << "}\n";

	//}

	return positions;
}




//since we fix a starting point but not an end point, the number of generated segments is equal to 
//the number of generated monomers
std::vector<std::vector<double>> random_walk(std::vector<double> &starting_end, int monomers, std::vector<std::vector<double>>& excluded_volumes)
{
	std::vector<std::vector<double>> generated_positions(monomers+1);//plus one because we're going to include
	// the starting point in the return vector of positions
	generated_positions[0] = starting_end;
	std::vector<double> old_position{ starting_end };
	std::vector<double> jump(3);
	std::vector<double> trial_position(3);

	for (int i{ 1 }; i <= monomers; i++) {
		sample_jump_direction(jump, 1);
		trial_position = vector_addition(old_position, jump);

		//generated_positions[i] = select_trial(trial_positions, excluded_volumes);
		generated_positions[i] = trial_position;

		excluded_volumes.push_back(generated_positions[i]);
		old_position = generated_positions[i];

	}
	return generated_positions;
}


std::vector<std::vector<double>> grow_chain(std::vector<double>& starting_end,
	std::vector<double>& ending_end, int l) {

	std::vector<std::vector<double>> positions(l + 1);
	std::vector<double> generated_position(3);
	std::vector<double> old_position{ starting_end };
	positions[0] = starting_end;
	int i_segments;


	//i is the monomer that is being grown. if you have 4 segments in total for example then you need
	//to add 3 monomers in addition to the endpoints so i=1,2,3. Note that including the endpoints we have
	//5 monomers in total in the flexible chain. 
	for (int i = 1; i <= l - 1; i++) {
		std::vector<double> trial_position(3);

		//std::cout <<"regrowth stage " << i << std::endl;
		if (l - i > 1) {
			i_segments = l - i;
			//if(i!=1)
			generated_position = rejection_sample(old_position, ending_end, i_segments);
			old_position = generated_position;
			positions[i] = generated_position;

		}
		else if (l - i == 1) {
			//std::cout << "Crankshaft regime" << std::endl;
			trial_position = crankshaft_insertion(positions[l - 2], ending_end);
			generated_position = trial_position;

			positions[l - 1] = generated_position;

		}
	}
	positions[l] = ending_end;

	return positions;
}
