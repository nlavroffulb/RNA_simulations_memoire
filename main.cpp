// the last dance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include "math_functions.h"
#include "helix_generation.h"
#include "simulations.h"
//#include "acceptance.h"
#include<time.h>
#include <fstream>
#include <chrono>
#include <algorithm>

//TO DO.  
////////////////////////////////////////////////////////////////////////////////////////

// appears to be a bug in which for a growth limit [a,b] between fixed ends, the weight of the (b-2)th monomer is often 0. this is the monomer right before crankshaft insertion. not sure what's going on but after 
// quite a bit of testing i can't find the bug, it seems to be working as intended. might come back to this later but for now i have to assume this is working as intended.
// 
// yamakawa function tests didn't work great. perhaps it's because im not testing in the long enough limit ? (number of 
// segments tends to infinity).
// 
// how does swivel work with zipped structures? before we assumed that the centre of mass was in the plane of 
// one of the base pairs however if we have an even number of base pairs e.g 4 this will not be the case, so is 
// it working correctly? probably not.
// 
// 
// can make metropolis structure more efficient. e,g only calculate search results, extensible structures when 
// a certain move has been accepted. after a series of rejections it should not be re calcualted.
// 
// might need to include normalization factor for yamakawa probability in the acceptance rules.
// 
// change zip so that we always sample the side of zip from which we can actually extend the helix. need to modify
// the acceptance rule potentially.
// 
// 
// acceptance rate is often nan. need to fix. 
// 
// 
// are we protected against structures which contain one of the terminal endpoints? need to check for all moves. 
// 
// unlink growth limits. sometimes we have growth limits that connect e.g [38,16] & [41,38]. should just be one growth limit.
// still need to include the extra probabilities which arise for each acceptance rule for example
// the probability of picking a zipped structure. 
// 
// structure id needs to be cleaned up i think. 
// // 
// when unlinking, is there any point in sampling alpha and beta?? the number of growth limits will be the same in either case.
// 
// limiting case for zip where we zip the two tail ends of the chain. am i protected against that?
// 
//  keep doing this thing for alpha  =  1 case where we reverse the order of the vector... is this necessary? when we do constraints
// we could just consider that w,x,y,z go from 3' to 5'. don't need to reverse it. remove the alpha dependency from w,x,y,z. then 
// there would be no need to reverse them.
// 
// 
////////////////////////////////////////////////////////////////////////////////////////

// 



// visualization of the structures that i make: different structures should have different colours.
// the monomers they enclose as well so we can easily follow the sequence.
// -- separate files for each structure. 
// -- different colours for connecting parts




int main()
{
    srand(time(NULL));
    std::int64_t init{ rand()};
    init = -1*init;
    ran2(&init);


    
    //polymer p0(70);
    double u1, u2, u3;
    int NMC{ 1000 },i{0};

    //hairpin_sim_varying_length(NMC, 110, true);
    link_unlink_hairpin_sim(NMC, 30, true);


     //   simple example of growing the polymer once and can output txt file to view it in ovito. can also measure
     //performance.
 

    //********************************************************************************
    //********************************************************************************
    //********************************************************************************

    //code to investigate e2e distance of middle monomer for polymer of length 71. 
    // results outputted to text file.
 
    //auto start = std::chrono::high_resolution_clock::now();
    //int N_sims{ 10000 };
    ////std::vector<double> sim_results(N_sims);
    //std::vector<double> sim_results(2*N_sims), temp(3);
    //std::vector<double> N_position{ 0,0,50 }, initial_position{ 0,0,0 };
    //int N_seg{ 120 };
    //int N_mid = N_seg / 2;
    //for (int i{ 0 }; i < N_sims; i++) {
    //    std::cout << i << std::endl;
    //    std::vector<std::vector<double>> p{ grow_chain(N_position,initial_position,N_seg) };
    //    temp = vector_subtraction(p[N_mid],initial_position);
    //    sim_results[2*i]=temp[2];
    //    temp = vector_subtraction(N_position, p[N_mid]);
    //    sim_results[2 * i + 1] = temp[2];

    //}

    //const char* path = "C:/Users/32471/Documents/Mémoire/testing yamakawa function/sim_results_x.txt";
    //std::ofstream file(path); //open in constructor

    //for (int k{ 0 }; k<sim_results.size(); k++) {
    //    file << sim_results[k] << std::endl;
    //}
    //std::ofstream vector("sim_results2.txt");
    //for (int k{ 0 }; k < N_sims; k++) {
    //    file << sim_results[k] << std::endl;
    //}

    //********************************************************************************
    //********************************************************************************
    //********************************************************************************



}
