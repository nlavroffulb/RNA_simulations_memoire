// the last dance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include "math_functions.h"
#include "helix_generation.h"
#include<time.h>
#include <fstream>
#include <chrono>
#include <algorithm>




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
    std::vector<std::vector<double>> p{ rosenbluth_random_walk(initial_position,50,ev0,energies,weights) };
    polymer p0(p);



    //METROPOLIS HASTINGS
    int NMC{ 1000 };

    int i{ 0 };
    double u1, u2;
    double link_acc{0.5}, unlink_acc{0.5};
    int link_counts{ 0 }, unlink_counts{ 0 };
    p0.link();
    while (i < NMC) {
        u1 = rand2(0, 1);
        u2 = rand2(0, 1);

        //LINK BRANCH
        if (u1 <= 0.5) {
            if (p0.linked()) { // if the system is already in the linked state we cannot re apply a link.
                link_counts++;
                i++;
                continue;
            }
            else {
                int alpha, beta, index;
                std::vector<int> s;
                p0.sample_link_region(s, alpha, beta, index);
                if (!p0.reject_link1(s, alpha, beta) && u2 <= link_acc) {
                    p0.link1(s, alpha, beta, index);
                    link_counts++;
                    std::cout << "link sucess" << std::endl;

                    //re calculate acceptance rules
                    link_acc = p0.acceptance_link0(1);
                }
                else {
                    unlink_counts++;
                }
            }            
        }

        //UNLINK BRANCH
        else {//u1 > 0.5 and less than 1.
            if (!p0.linked()) {// if the system is already in the unlinked state we cannot re apply an unlink move.
                unlink_counts++;
                i++;
                continue;
            }
            else {
                if (u2 <= unlink_acc) {
                    p0.unlink();// attempt an unlink move. it is never rejected but we still check.
                    unlink_counts++;
                    std::cout << "unlink sucess" << std::endl;

                    //recalculate acceptance probabilities
                    unlink_acc = p0.acceptance_link0(0);
                }
                else {
                    link_counts++;
                }
            }
        }
        i++;
    }
    std::cout << "Unlink state: " << unlink_counts << std::endl;
    std::cout << "Link state: " << link_counts << std::endl;
    
    double time{ 0.0 };
    int trials{ 1000 };
    int alpha, beta, s_index;
    std::vector<int> s;
    int l{ 0 };
    for (int i{ 0 }; i < trials; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        int c{ 0 };

        std::cout << i << std::endl;
        polymer p0(p);

        std::cout << "1st link move " << std::endl;
        p0.sample_link_region(s, alpha, beta,s_index);
        if (!p0.reject_link1(s, alpha, beta)) {
            p0.link1(s, alpha, beta,s_index);
            std::cout << "SUCCESS" << std::endl;
        }
        else {
            std::cout << "link rejected" << std::endl;
        }
        p0.output_for_ovito4("step1");

        std::cout << "2nd link move " << std::endl;
        p0.sample_link_region(s, alpha, beta, s_index);
        if (!p0.reject_link1(s, alpha, beta)) {
            p0.link1(s, alpha, beta, s_index);
            std::cout << "SUCCESS" << std::endl;



        }
        else {
            std::cout << "link rejected" << std::endl;
        }

        p0.output_for_ovito4("step2");
        //p0.unlink();
        std::cout << "3rd link move " << std::endl;
        p0.sample_link_region(s, alpha, beta, s_index);
        if (!p0.reject_link1(s, alpha, beta)) {
            p0.link1(s, alpha, beta, s_index);
            std::cout << "SUCCESS" << std::endl;



        }
        else {
            std::cout << "link rejected" << std::endl;
        }
        p0.output_for_ovito4("step3");

        std::cout << "zip move" << std::endl;
        bool move_success;
        p0.zip(move_success);
        p0.output_for_ovito4("step4");

        if (move_success) {
            c++;
        }

        std::cout << "unzip move" << std::endl;
        p0.unzip(move_success);
        if (move_success) {
            c++;
        }
        p0.output_for_ovito4("step5");

        std::cout << "2nd zip move" << std::endl;
        p0.zip(move_success);
        if (move_success) {
            c++;
        }
        p0.output_for_ovito4("step6");
        if (c == 3) {
            std::cout << "stop";
        }
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


    //std::vector<std::vector<double>> s_x, s_y,s_x1,s_y1;
    //std::vector<double> u{1,0,0}, v{0,1,0}, u1{-1,0,0}, v1{0,-1,0};
    //std::vector<double> op{ -0.38716,2,0.0716736 };
    //sample_helix_vectors(u, v);
    //for (auto i : u) {
    //    std::cout << i << std::endl;
    //}
    //for (auto i : v) {
    //    std::cout << i << std::endl;
    //}
    //generate(3, v, u, initial_position, s_x, s_y);
    //print_2d_doub_vec(s_x);
    //print_2d_doub_vec(s_y);

    //generate(3, v1, u1, op, s_x1, s_y1);
    //print_2d_doub_vec(s_x1);
    //print_2d_doub_vec(s_y1);

    //p0.link();
    //p0.output_for_ovito4("smow0");

    //helix_struct* my_struct{ p0.sample_double_helix() };
    //p0.add_bp_to_helix(my_struct, 1);
    //p0.output_for_ovito4("smow1");


    //for (int i{ 0 }; i < 3; i++) {
    //    print_1d_doub_vec
    //}
    //for (auto i : u) {
    //    std::cout << i << std::endl;
    //}
    //for (auto i : v) {
    //    std::cout << i << std::endl;
    //}


}
