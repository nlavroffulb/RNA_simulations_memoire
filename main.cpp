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


    //std::cout << partition_ni_unlinked(3) << std::endl;

    //std::vector<double> u{ 1,0,0 }, v{ 0,1,0 }, o{ 0,0,0 };
    //sample_helix_vectors(u, v);
    //int n{ 10 };
    //std::vector<std::vector<double>> s_x0(n), s_y0(n), s_x(n), s_y(n);

    //generate(n, v, u, o, s_x0, s_y0);
    //std::cout << std::endl;
    //print_2d_doub_vec(s_x0);
    //print_2d_doub_vec(s_y0);
    //std::cout << dist_2_points3d(s_y0[0], s_y0[1])<<std::endl;
    //std::cout << consecutive_mon_sep() << std::endl;
    //std::cout << max_ring_sep() << std::endl;
    //std::cout << vector_modulus(v) << std::endl;

    //std::cout << "reverse generate " << std::endl;
    //reverse_generate(n, v, u, s_x0.back(), s_y0.front(), s_x, s_y);
    //std::cout << std::endl;
    //print_2d_doub_vec(s_x);
    //print_2d_doub_vec(s_y);


    ////METROPOLIS HASTINGS
    //int NMC{ 1000 };

    //int i{ 0 };
    //double u1, u2;
    //double link_acc{0.5}, unlink_acc{0.5};
    //int link_counts{ 0 }, unlink_counts{ 0 };
    //p0.link();
    //while (i < NMC) {
    //    u1 = rand2(0, 1);
    //    u2 = rand2(0, 1);

    //    //LINK BRANCH
    //    if (u1 <= 0.5) {
    //        if (p0.linked()) { // if the system is already in the linked state we cannot re apply a link.
    //            link_counts++;
    //            i++;
    //            continue;
    //        }
    //        else {
    //            int alpha, beta;
    //            std::vector<int> s;
    //            p0.sample_link_region(s, alpha, beta);
    //            if (!p0.reject_link1(s, alpha, beta) && u2 <= link_acc) {
    //                p0.link1(s, alpha, beta);
    //                link_counts++;
    //                std::cout << "link sucess" << std::endl;

    //                //re calculate acceptance rules
    //                link_acc = p0.acceptance_link0(1);
    //            }
    //            else {
    //                unlink_counts++;
    //            }
    //        }            
    //    }

    //    //UNLINK BRANCH
    //    else {//u1 > 0.5 and less than 1.
    //        if (!p0.linked()) {// if the system is already in the unlinked state we cannot re apply an unlink move.
    //            unlink_counts++;
    //            i++;
    //            continue;
    //        }
    //        else {
    //            if (u2 <= unlink_acc) {
    //                p0.unlink();// attempt an unlink move. it is never rejected but we still check.
    //                unlink_counts++;
    //                std::cout << "unlink sucess" << std::endl;

    //                //recalculate acceptance probabilities
    //                unlink_acc = p0.acceptance_link0(0);
    //            }
    //            else {
    //                link_counts++;
    //            }
    //        }
    //    }
    //    i++;
    //}
    //std::cout << "Unlink state: " << unlink_counts << std::endl;
    //std::cout << "Link state: " << link_counts << std::endl;
    
    double time{ 0.0 };
    int trials{ 10 };
    int alpha, beta, s_index;
    std::vector<int> s;
    int l{ 0 };
    for (int i{ 0 }; i < trials; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        int c{ 0 };
        bool success{ 0 };
        std::cout << i << std::endl;
        polymer p0(50);

        std::cout << "1st link move " << std::endl;
        p0.sample_link_region(s, alpha, beta,s_index);
        if (!p0.reject_link1(s, alpha, beta)) {
            p0.link(s, alpha, beta,s_index);
            std::cout << "SUCCESS" << std::endl;
            p0.update_positions();

        }
        else {
            std::cout << "link rejected" << std::endl;
        }
        p0.output_for_ovito4("step1");

        std::cout << "2nd link move " << std::endl;
        p0.sample_link_region(s, alpha, beta, s_index);
        if (!p0.reject_link1(s, alpha, beta)) {
            p0.link(s, alpha, beta, s_index);
            std::cout << "SUCCESS" << std::endl;

            p0.update_positions();

        }
        else {
            std::cout << "link rejected" << std::endl;
        }

        p0.output_for_ovito4("step2");
        //p0.unlink();
        std::cout << "3rd link move " << std::endl;
        p0.sample_link_region(s, alpha, beta, s_index);
        if (!p0.reject_link1(s, alpha, beta)) {
            p0.link(s, alpha, beta, s_index);
            std::cout << "SUCCESS" << std::endl;
            p0.update_positions();



        }
        else {
            std::cout << "link rejected" << std::endl;
        }
        p0.output_for_ovito4("step3");
        
        p0.swivel1(success);
        p0.output_for_ovito4("step4");
        //std::cout << "zip move" << std::endl;
        //bool move_success;
        //p0.zip(move_success);
        //p0.output_for_ovito4("step4");

        //if (move_success) {
        //    p0.update_positions();

        //    c++;
        //}

        //std::cout << "unzip move" << std::endl;
        //p0.unzip(move_success);
        //if (move_success) {
        //    p0.update_positions();

        //    c++;
        //}
        //p0.output_for_ovito4("step5");

        //std::cout << "2nd zip move" << std::endl;
        //p0.zip(move_success);
        //if (move_success) {
        //    p0.update_positions();

        //    c++;
        //}
        //p0.output_for_ovito4("step6");
        //if (c == 3) {
        //    std::cout << "stop";
        //}
        //p0.unzip();



        //std::cout << "3rd link move " << std::endl;
        //p0.generalized_link();
        //std::cout << "UNLINK" << std::endl;
        //p0.unlink();
        //p0.unlink();
        //p0.generalized_link();
        //std::cout << "4th link move " << std::endl;

        //p0.generalized_link();
        //p0.output_for_ovito4();
        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << duration.count() << std::endl;
        time += duration.count();

    }
    std::cout << "3 successful links " << l << std::endl;
    std::cout << "average time " << time / static_cast<double>(trials) << std::endl;



}
