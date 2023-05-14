//#include "acceptance.h"
//#include "polymer_generation.h"
////#include "unbound_section_class.h"
//const int NMC{ 10000 };
//int trials{ 1000 };
//double partition_ni_unlinked(int length) {
//
//	int i{ 0 };
//	std::vector<double> old_position{ 0,0,0 };
//	std::vector<double> jump(3);
//	std::vector<double> weights(length);
//	std::vector<double> energies(length);
//	std::vector<double> generated_position(3);
//
//	std::vector <std::vector<double>> excluded_volume(length);
//
//	std::vector<std::vector<double>> trial_positions(trials);
//	std::vector<double> trial_position(3);
//	std::vector<double> NMC_weights(NMC);
//	while (i < NMC) {
//
//		for (int m{ 0 }; m < length;m++) {
//			if (m == 0) {
//				sample_jump_direction(jump, 1);
//				old_position = vector_addition(old_position, jump);
//				energies[0] = 0;
//				weights[0] = trials * exp(-energies[0]);
//				excluded_volume[0] = old_position;
//				
//			}
//			else {
//				for (int j{ 0 }; j < trials; j++) {
//					sample_jump_direction(jump, 1);
//					trial_positions[j] = vector_addition(old_position, jump);;
//				}
//				generated_position = select_trial(trial_positions, excluded_volume, weights[m], energies[m]);
//				old_position = generated_position;
//				excluded_volume[m] = old_position;
//			}
//
//		}
//		NMC_weights[i] = rosenbluth_factor(weights, trials)/NMC;
//
//		i++;
//		std::cout << i << std::endl;
//	}
//	return sum_of_elements(NMC_weights);
//
//}
//
//double partition_ni_linked(int length)
//{
//
//	return 0.0;
//}
//
