#include "simulations.h"
#include "math_functions.h"

// this is the simulation where only one link move is allowed at a time. corresponds to the formation of a hairpin loop.
void hairpin_sim(polymer* p0, int NMC)
{
    //polymer p0(N_monomers, rosenbluth);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    while (i < NMC) {
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            p0->link_move();
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            p0->unlink_move();
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }
        i++;
    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

}

void determine_g0(polymer* p0, int NMC)
{
    //polymer p0(N_monomers, true);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{0}, uc_0{0};

    double init_prefactor{ 0.05 };
    p0->set_link_acceptance_prefactor(init_prefactor);
    double r;
    std::vector<double> prefactors;
    while (i < NMC) {

        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            std::cout << "link branch" << std::endl;
            p0->link_move();
        }
        else {
            std::cout << "unlink branch" << std::endl;

            p0->unlink_move();
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i+1) % 5  == 0 && i!=0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;
            p0->set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0->get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;

        }
        i++;
    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

    print_1d_doub_vec(prefactors);

}

void linker_weight_sim(polymer* p0, int NMC)
{

    //polymer p0(N_monomers, true);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    std::vector<double> linker_weights;
    while (i < NMC) {
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            std::cout << "link branch" << std::endl;
            p0->link_move();
        }
        else {
            std::cout << "unlink branch" << std::endl;

            p0->unlink_move();
        }

        if (p0->linked() == true) {
            bound_counts++;
            linker_weights.push_back(p0->get_hairpin_weight());
        }
        else {
            unbound_counts++;
        }
        i++;
    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;
    
}

void varying_helix_length_sim(polymer* p0, int NMC)
{
    int i{ 0 },j{0};
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    int helix_max_length{ 10 };
    int N{ p0->get_num_monomers() };

    std::vector<double> link_probabilities(helix_max_length - 3);
    double link_prob;
    //p0->set_link_acceptance_prefactor(1);

    while (i < helix_max_length-3) {
        
        p0->set_search_results({ 0,2 + i,N - 3 - i,N - 1 });
        p0->force_unbound_state();//always begin in unlinked state
        j = 0;
        bound_counts = 0;
        unbound_counts = 0;
        if (i == 4) {
            std::cout << "stop" << std::endl;
        }
        while (j < NMC) {
            std::cout << "i " << i << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0->link_move();
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0->unlink_move();
            }

            if (p0->linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            j++;
        }
        if (unbound_counts == 0) {
            link_prob = 1;
        }
        else {
            link_prob = static_cast<double>(bound_counts) / static_cast<double>(bound_counts + unbound_counts);
        }
        link_probabilities[helix_max_length - 4 - i] = link_prob;

        i++;

    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;
    print_1d_doub_vec(link_probabilities);
}

void vary_dangle_sim(polymer* p0, int NMC)
{

    int max_dangle_length{ 10 };
    std::vector<double> helix_weights(max_dangle_length), helix_weights_l(NMC);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int N{ p0->get_num_monomers() };

    for (int l{ 0 }; l < max_dangle_length; l++) {
        p0->set_search_results({ 0 + l,2 + l,N - 3,N - 1 });
        p0->force_unbound_state();//always begin in unlinked state
        i = 0;
        bound_counts = 0;
        unbound_counts = 0;
        while (i < NMC) {
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0->link_move();
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0->unlink_move();
            }

            if (p0->linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            helix_weights_l[i] = p0->get_helix_weight();
            i++;
        }
        helix_weights[l] = sum_of_elements(helix_weights_l) / NMC;

    }
    print_1d_doub_vec(helix_weights);
    //std::cout << "Unlink state: " << unbound_counts << std::endl;
    //std::cout << "Link state: " << bound_counts << std::endl;

}

void kissing_hairpin_swivel(polymer* p0, int NMC, bool swivel_allowed)
{
    double swivel_switch;
    swivel_allowed == true ? swivel_switch = 2.0/3.0 : swivel_switch = 0.5;
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    int kissing_hairpin_count{ 0 };

    while (i < NMC) {
        u_branch = rand2(0, swivel_switch);
        if (u_branch < 1.0/3.0) {// link_branch
            //std::cout << "link branch" << std::endl;
            p0->link_move();
        }
        else if(u_branch <2.0/3.0) {
            //std::cout << "unlink branch" << std::endl;

            p0->unlink_move();
        }
        else if (u_branch < 1) {
            p0->swivel_move();
        }
        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }
        if (p0->get_num_helices() == 3) {
            kissing_hairpin_count++;
        }
        i++;
    }

    double P_kissing_hairpin{ static_cast<double>(kissing_hairpin_count) / static_cast<double>(NMC) };
}

