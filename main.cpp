// the last dance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include "math_functions.h"
#include "helix_generation.h"
//#include "acceptance.h"
#include<time.h>
#include <fstream>
#include <chrono>
#include <algorithm>

//#include "ran2.h"

// average time for 2 consecutive link moves: 1937 ms (almost 2 seconds yikes). 4/4 15:00. need to get this down a lot.
// 1170ms 5/4 20:00 KEEP GOING !!!
// 830ms 5/4 20:30 RAHHHHHHHHHH!!!

// average time for 3 consecutive link moves 27/4: 700ms. quite good.

// time for 2 links and 1 unlink
// 1356ms 5/4 21:00
// 
// 
// 
// 
//TO DO.  
////////////////////////////////////////////////////////////////////////////////////////

// global constants header file.
// 
// acceptance rate is often nan. need to fix. 
// 
// do we always reject links if there are no structures left to try?
// 
// are we protected against structures which contain one of the terminal endpoints? need to check for all moves. 
// 
// unlink growth limits. sometimes we have growth limits that connect e.g [38,16] & [41,38]. should just be one growth limit.
// still need to include the extra probabilities which arise for each acceptance rule for example
// the probability of picking a zipped structure. 
// 
// structure id needs to be cleaned up i think. 
// 
// when writing acceptance rules, will be useful to use growth limits functions? 
// 
// when unlinking, is there any point in sampling alpha and beta?? the number of growth limits will be the same in either case.
// 
// limiting case for zip where we zip the two tail ends of the chain. am i protected against that?
// 
//  keep doing this thing for alpha  =  1 case where we reverse the order of the vector... is this necessary? when we do constraints
// we could just consider that w,x,y,z go from 3' to 5'. don't need to reverse it. remove the alpha dependency from w,x,y,z. then 
// there would be no need to reverse them.
// 
// for metropolis hastings: 
// ------------all our MC move functions shouldn't directly change the positions of the monomers because even if a move is successful
// it could be rejected because of the acceptance probability.
// 
// 
// build swivel
// 
// work on decoupling the various parts of each move to accommodate physical constraints and the acceptance probability.
// 
// for our simple link unlink simulation, could initially store the weights for the whole chain then modify as we go along (regrowing certain sections
// by random walk). create different objects for fixed end 2 end growth sections. for this case its simple, we'll only ever have 1 fixed end to 
// end growth section, but good to see how this will scale for more complicated situations.
// 
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
    //double jump{ rand2(0,1) };
    std::vector<double> initial_position{0,0,0};
    std::vector<double> N_position{ 0,0,30 }, weights, energies;
    std::vector<std::vector<double>> ev0;
    //std::vector<std::vector<double>> p{ rosenbluth_random_walk(initial_position,50,ev0,energies,weights) };
    //polymer p0(50);


    

    polymer p0(70);
    double u1, u2, u3;
    int NMC{ 200 },i{0};

    while (i<NMC){
        p0.neighbouring_linkers();
        //std::cout << "Round " << i << std::endl;
        u1 = rand2(0, 1);// select link/unlink, zip/unzip or swivel branches
        u2 = rand2(0, 1);// select forward or backwards move
        u3 = rand2(0, 1);// reject based on acceptance probability.
        //u1 = rand2(0, 1)*0.66;;
        //u2 = 0.5;
        bool success;

        std::string filename{ "step" + std::to_string(i) };
        p0.output_for_ovito4(filename);

        if (u1 < 0.3333)//link unlink branch
        {
            int alpha, beta, struct_index;
            std::vector<int> s;
            if (u2 <= 0.5) {
                //link subbranch
                std::cout << "link branch" << std::endl;

                p0.sample_link_region(s, alpha, beta, struct_index);
                if (!p0.reject_link(s, alpha, beta)) {//physical constraints
                    //p0.print_monomer_positions();
                    p0.link(s, alpha, beta, struct_index);

                    if (u3 < p0.link_acceptance(1)) {// acceptance probability
                        print_1d_int_vec(s);
                        std::cout << "LINK ACCEPTED" << std::endl;
                        p0.link_update(struct_index);

                        print_1d_int_vec(s);
                    }
                }


            }
            else {
                //unlink subbranch
                std::cout << "unlink branch" << std::endl;

                if (!p0.linked()) {// if there are no links then we cannot unlink.
                    i++;
                    continue;
                }
                int s_index;
                p0.sample_unlink_region(s_index);
                if (s_index == -1) {// we sampled a zipped structure. we cannot unlink zipped structures.
                    i++;
                    continue;
                }
                else {
                    p0.unlink(s_index);
                    if (u3 < p0.link_acceptance(0)) {
                        p0.unlink_update(s_index);
                        std::cout << "UNLINK ACCEPTED" << std::endl;

                    }
                }


            }

        }
        else if (u1 < 0.6666) { // zip unzip branch
            int side, s_index;
            if (u2 < 0.5) {//zip branch
                std::cout << "zip branch" << std::endl;
                if (p0.linked()) {
                    p0.zip(success, side, s_index);
                }
                else {
                    p0.zip(success, side, s_index);

                }
                if (success == 1) { // success is dependent on whether there are any suitable double helixes which already exist and can be extended (ie zipped).
                    if (u3 < p0.zip_acceptance(1,side)) {
                        p0.zip_update(s_index,side);
                        std::cout << "ZIP ACCEPTED" << std::endl;

                    }
                }
            }
            else {//unzip branch
                std::cout << "unzip branch" << std::endl;

                p0.unzip(success, side,s_index);
                if (success == 1) {
                    if (u3 < p0.zip_acceptance(0, side)) {
                        p0.unzip_update(s_index,side);
                        std::cout << "UNZIP ACCEPTED" << std::endl;

                    }
                }
            }

        }
        else { // swivel branch
            helix_struct* dh{ p0.sample_double_helix(success) };
            if (success == 0) {// there were no existing double helix structures.
                i++;
                continue;
            }
            double tau{ 0.5 };
            double theta{ rand2(0,tau * atan(1) * 2) };
            int kappa{ static_cast<int>(rand2(0,2)) };
            theta = theta * pow(-1, kappa);

            if (u2 < 0.3333) {// spin branch. corresponds to changing the v vector of the helix
                //std::cout << "v swivel branch" << std::endl;

                if (p0.reject_spin(dh, theta)) {
                    //std::cout << "v swivel rejected" << std::endl;
                    continue;
                }
                else {
                    std::vector<double> rotation(3);
                    p0.spin(dh, theta, rotation);

                    if (u3 < p0.swivel_acceptance()) {
                        p0.spin_update(dh,rotation,theta);
                        std::cout << "V SWIVEL ACCEPTED" << std::endl;

                    }
                }

            }
            else if (u2 < 0.6666) {// corkscrew branch. corresponds to changing the u vector of the helix
                //std::cout << "u swivel branch" << std::endl;

                if (p0.reject_corkscrew(dh, theta)) {
                    //std::cout << "u swivel rejected" << std::endl;
                    continue;
                }
                else {
                    p0.corkscrew(dh, theta);

                    if (u3 < p0.swivel_acceptance()) {
                        p0.corkscrew_update(dh);
                        std::cout << "U SWIVEL ACCEPTED" << std::endl;

                    }
                }
            }
            else {// translation branch. displace the helix.
                //std::cout << "origin translation branch" << std::endl;

                double translation_distance{ 0.5 };
                std::vector<double> translation(3);
                sample_jump_direction(translation, translation_distance);

                if (p0.reject_translate(dh, translation)) {
                    //std::cout << "origin translation rejected" << std::endl;
                    continue;
                }
                else {
                    p0.translate(dh, translation);

                    if (u3 < p0.swivel_acceptance()) {
                        p0.translate_update(dh, translation);
                        std::cout << "O TRANSLATION ACCEPTED" << std::endl;

                    }
                }
            }
        }
        i++;
    }




}
