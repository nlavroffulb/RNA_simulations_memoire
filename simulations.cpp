#include "simulations.h"
#include "math_functions.h"

// this is the simulation where only one link move is allowed at a time. corresponds to the formation of a hairpin loop.
void hairpin_sim(int NMC, int N_monomers, bool rosenbluth)
{
    polymer p0(N_monomers, rosenbluth);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    while (i < NMC) {
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            std::cout << "link branch" << std::endl;
            p0.link_move();
        }
        else {
            std::cout << "unlink branch" << std::endl;

            p0.unlink_move();
        }

        if (p0.linked() == true) {
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

        if ((i+1) % 20  == 0 && i!=0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;
            p0->set_link_acceptance_prefactor(0, r);
            bc_0 = bound_counts, uc_0 = unbound_counts;
            prefactors.push_back(p0->get_link_acceptance_prefactor());
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
        p0->force_unlink_move();//always begin in unlinked state
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

