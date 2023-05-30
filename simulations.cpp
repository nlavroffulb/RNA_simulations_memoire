#include "simulations.h"
#include "polymer_class.h"
#include "math_functions.h"

// this is the simulation where only one link move is allowed at a time. corresponds to the formation of a hairpin loop.
void link_unlink_hairpin_sim(int NMC, int N_monomers, bool rosenbluth)
{
    //int N_monomers{ 50 };
    polymer p0(N_monomers, rosenbluth);
    //METROPOLIS HASTINGS
    //int NMC{ 1000 };

    int i{ 0 };
    double u1, u2, u3;
    double link_acc{ 0.5 }, unlink_acc{ 0.5 };
    int link_counts{ 0 }, unlink_counts{ 0 };

    int alpha, beta, struct_index;

    while (i < NMC) {
        u1 = rand2(0, 1);
        u2 = rand2(0, 1);
        u3 = rand2(0, 1);

        //LINK BRANCH
        if (u1 <= 0.5) {
            int alpha, beta, struct_index;
            std::vector<int> s;

            if (p0.linked()) { // if the system is already in the linked state we cannot re apply a link.
                link_counts++;
                i++;
                continue;
            }
            else {
                p0.sample_link_region(s, alpha, beta, struct_index);
                if (!p0.reject_link(s, alpha, beta) && u2 <= link_acc) {
                    p0.link(s, alpha, beta, struct_index);
                    std::cout << "link sucess" << std::endl;

                    if (u3 < p0.link_acceptance(1)) {// acceptance probability
                        print_1d_int_vec(s);
                        std::cout << "LINK ACCEPTED" << std::endl;
                        p0.link_update(struct_index);
                        link_counts++;
                        //print_1d_int_vec(s);
                    }
                    else {
                        unlink_counts++;
                    }
                }
                else {
                    unlink_counts++;
                }
            }
        }

        //UNLINK BRANCH
        else {//u1 > 0.5 and less than 1.
            std::cout << "unlink branch" << std::endl;

            if (!p0.linked()) {// if there are no links then we cannot unlink.
                i++;
                unlink_counts++;
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
                    unlink_counts++;
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

}



//while (i<NMC){
//    std::cout << "frame " << i << std::endl;
//    p0.neighbouring_linkers();

//    //std::cout << "Round " << i << std::endl;
//    u1 = rand2(0, 1);// select link/unlink, zip/unzip or swivel branches
//    u2 = rand2(0, 1);// select forward or backwards move
//    u3 = rand2(0, 1);// reject based on acceptance probability.
//    //u1 = rand2(0, 1)*0.66;;
//    //u2 = 0.5;
//    bool success;

//    std::string filename{ "step" + std::to_string(i) };
//    //p0.output_for_ovito4(filename);
//    p0.ovito_bonds(filename);
//    if (u1 < 0.3333)//link unlink branch
//    {
//        int alpha, beta, struct_index;
//        std::vector<int> s;
//        if (u2 <= 0.5) {
//            //link subbranch
//            std::cout << "link branch" << std::endl;

//            p0.sample_link_region(s, alpha, beta, struct_index);
//            if (!p0.reject_link(s, alpha, beta)) {//physical constraints
//                //p0.print_monomer_positions();
//                p0.link(s, alpha, beta, struct_index);

//                if (u3 < p0.link_acceptance(1)) {// acceptance probability
//                    print_1d_int_vec(s);
//                    std::cout << "LINK ACCEPTED" << std::endl;
//                    p0.link_update(struct_index);

//                    //print_1d_int_vec(s);
//                }
//            }


//        }
//        else {
//            //unlink subbranch
//            std::cout << "unlink branch" << std::endl;

//            if (!p0.linked()) {// if there are no links then we cannot unlink.
//                i++;
//                continue;
//            }
//            int s_index;
//            p0.sample_unlink_region(s_index);
//            if (s_index == -1) {// we sampled a zipped structure. we cannot unlink zipped structures.
//                i++;
//                continue;
//            }
//            else {
//                p0.unlink(s_index);
//                if (u3 < p0.link_acceptance(0)) {
//                    p0.unlink_update(s_index);
//                    std::cout << "UNLINK ACCEPTED" << std::endl;

//                }
//            }


//        }

//    }
//    else if (u1 < 0.6666) { // zip unzip branch
//        int side, s_index;
//        if (u2 < 0.5) {//zip branch
//            std::cout << "zip branch" << std::endl;
//            if (p0.linked()) {
//                p0.zip(success, side, s_index);
//            }
//            else {
//                p0.zip(success, side, s_index);

//            }
//            if (success == 1) { // success is dependent on whether there are any suitable double helixes which already exist and can be extended (ie zipped).
//                if (u3 < p0.zip_acceptance(1,side)) {
//                    p0.zip_update(s_index,side);
//                    std::cout << "ZIP ACCEPTED" << std::endl;

//                }
//            }
//        }
//        else {//unzip branch
//            std::cout << "unzip branch" << std::endl;

//            p0.unzip(success, side,s_index);
//            if (success == 1) {
//                if (u3 < p0.zip_acceptance(0, side)) {
//                    p0.unzip_update(s_index,side);
//                    std::cout << "UNZIP ACCEPTED" << std::endl;

//                }
//            }
//        }

//    }
//    else { // swivel branch
//        helix_struct* dh{ p0.sample_double_helix(success) };
//        if (success == 0) {// there were no existing double helix structures.
//            i++;
//            continue;
//        }
//        double tau{ 0.5 };
//        double theta{ rand2(0,tau * atan(1) * 2) };
//        int kappa{ static_cast<int>(rand2(0,2)) };
//        theta = theta * pow(-1, kappa);

//        if (u2 < 0.3333) {// spin branch. corresponds to changing the v vector of the helix
//            //std::cout << "v swivel branch" << std::endl;

//            if (p0.reject_spin(dh, theta)) {
//                //std::cout << "v swivel rejected" << std::endl;
//                continue;
//            }
//            else {
//                std::vector<double> rotation(3);
//                p0.spin(dh, theta, rotation);

//                if (u3 < p0.swivel_acceptance()) {
//                    p0.spin_update(dh,rotation,theta);
//                    std::cout << "V SWIVEL ACCEPTED" << std::endl;

//                }
//            }

//        }
//        else if (u2 < 0.6666) {// corkscrew branch. corresponds to changing the u vector of the helix
//            //std::cout << "u swivel branch" << std::endl;

//            if (p0.reject_corkscrew(dh, theta)) {
//                //std::cout << "u swivel rejected" << std::endl;
//                continue;
//            }
//            else {
//                p0.corkscrew(dh, theta);

//                if (u3 < p0.swivel_acceptance()) {
//                    p0.corkscrew_update(dh);
//                    std::cout << "U SWIVEL ACCEPTED" << std::endl;

//                }
//            }
//        }
//        else {// translation branch. displace the helix.
//            //std::cout << "origin translation branch" << std::endl;

//            double translation_distance{ 0.5 };
//            std::vector<double> translation(3);
//            sample_jump_direction(translation, translation_distance);

//            if (p0.reject_translate(dh, translation)) {
//                //std::cout << "origin translation rejected" << std::endl;
//                continue;
//            }
//            else {
//                p0.translate(dh, translation);

//                if (u3 < p0.swivel_acceptance()) {
//                    p0.translate_update(dh, translation);
//                    std::cout << "O TRANSLATION ACCEPTED" << std::endl;

//                }
//            }
//        }
//    }
//    i++;
//}

