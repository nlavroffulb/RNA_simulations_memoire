#include "simulations.h"
#include "math_functions.h"

// this is the simulation where only one link move is allowed at a time. corresponds to the formation of a hairpin loop.
void hairpin_sim(polymer* p0, int NMC)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    
    while (i < NMC) {
        u_branch = rand2(0, 1);
        //if (p0->get_num_helices() == 2) {
        //    std::cout << "stop " << std::endl;
        //}
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0->link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0->unlink_move(m);
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

// might need a new strategy for convergence. it does seem correct, just the convergence rate isn't good.
void determine_prefactor(polymer* p0, int NMC, int update_rate, std::string filename)
{
    //polymer p0(N_monomers, true);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{0}, uc_0{0};

    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };
    p0->set_link_acceptance_prefactor(init_prefactor);
    double r{0};
    std::vector<double> prefactors;
    prefactors.push_back(init_prefactor);
    int cc{ 0 };
    while (i < NMC) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    prefactors.push_back(0.000001);
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0->link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0->unlink_move(m);
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i+1) % update_rate  == 0 && i!=0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;
          
            p0->set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0->get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0->get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

    double_vector_to_txt("Results/determine_g0.txt",prefactors);



}

double get_prefactor(std::vector<int> structure, int N_monomers, int NMC, int update_rate, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };

    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };

    polymer p0(N_monomers, true);
    p0.set_search_results(structure);
    p0.set_link_acceptance_prefactor(init_prefactor);
    double r{ 0 };
    std::vector<double> prefactors;
    prefactors.push_back(init_prefactor);
    int cc{ 0 };

    int t{ 0 };
    while (i < NMC) {
        if ((i + 1) % update_rate == 0 && i != 0) {
            t = i + 1;
        }
        // convergence condition.
        if (i==t && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    //prefactors.push_back(0.000001);
                    return prefactors.back();
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i + 1) % update_rate == 0 && i != 0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;

            p0.set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0.get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0.get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }

    return prefactors.back();
}

double correlation_length(double acceptance_prefactor, int N_monomers, std::vector<int> structure, int NMC, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };    


    std::vector<int> correlation_lengths;
    bool link_state0{ false }, link_state{false};
    int last_switch{0};


    // a run to calculate the correlation time for a given prefactor. this is the average time that
    // the system takes to switch from the linked to the unlinked configuration.
    polymer p0(N_monomers, true);
    p0.set_link_acceptance_prefactor(acceptance_prefactor);
    p0.set_search_results(structure);
    while (i < NMC) {
        std::cout << "iteration " << i << " out of " << NMC << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            if (link_state0 == false) {
                link_state0 = true;
                correlation_lengths.push_back(i - last_switch+1);
                last_switch = i+1;
            }
            bound_counts++;
        }
        else {
            if (link_state0 == true) {
                link_state0 = false;
                correlation_lengths.push_back(i - last_switch+1);
                last_switch = i+1;
            }
            unbound_counts++;
        }
        i++;
    }

    // calculate correlation length
    int_vector_to_txt("Results/correlation_length" + filename + ".txt", correlation_lengths);

    double avg_correlation{ static_cast<double>(sum_of_elements(correlation_lengths)) / static_cast<double>(correlation_lengths.size()) };
    int update_rate{ static_cast<int>(std::round(avg_correlation)) };

    // we want to make sure that we do a number of trials that is at least 10 times the correlation length 
    // so that there will be at least 10 updates of the prefactor.
    int NMC2{ NMC };
    if (avg_correlation / static_cast<int>(NMC) < 10) {
        NMC2 = 10 * NMC2;
    }
    int bc_0{ 0 }, uc_0{ 0 };
    double epsilon{ 2 };// for convergence
    int cc{ 0 };
    double r{ 0 };
    std::vector<double> prefactors{ acceptance_prefactor };

    i = 0;
    // another run to get a more precise prefactor. compared to the previous run, the re calculation of the 
    // prefactor happens at intervals equal to the calculated correlation time.
    while (i < NMC2) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    //prefactors.push_back(0.000001);
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;
            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i + 1) % update_rate == 0 && i != 0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;

            p0.set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0.get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0.get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }

    return prefactors.back();
    double_vector_to_txt("Results/g0_correlation_length" + filename + ".txt", prefactors);

}

void linker_weight_sim(int N_monomers, int NMC, std::string path)
{

    int max_dangle_length{ 20 };

    //data 
    std::vector<double> link_probability(max_dangle_length), link_acceptances, unlink_acceptances, prefactors(max_dangle_length);

    std::vector<double> linker_weights(max_dangle_length), helix_weights(max_dangle_length),
        linker_weights_l, helix_weights_l, physical_rejections(max_dangle_length),
        old_linker_weights, new_linker_weights, old_dangle_weights, new_dangle_weights;
    double overstretch_rejections{ 0 };

    bool skip_trial{ false };


    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    std::vector<double> old_weights, new_weights;
    //p0->force_unbound_state();

    int N0{ N_monomers }, N{ N0 };
    std::vector<int> helix(4);

    for (int l{ 0 }; l < max_dangle_length; l++) {

        // the tail is going to get longer such that the total number of monomers in the polymer increases each time.
        N = N0 + 2 * l;
        polymer p0(N, true);
        helix = { 0 + l, 2 + l, N - 3 - l, N - 1 - l };
        p0.set_search_results(helix);



        prefactors[l] = get_prefactor(helix, N, NMC, 100);
        prefactors[l] = correlation_length(prefactors[l], N, helix, NMC, std::to_string(l));
        p0.set_link_acceptance_prefactor(prefactors[l]);

        std::cout << "dangle length " << l << std::endl;

        i = 0;
        bound_counts = 0;
        unbound_counts = 0;

        helix_weights_l.clear();
        linker_weights_l.clear();
        old_weights.clear();
        new_weights.clear();
        link_acceptances.clear();
        unlink_acceptances.clear();
        overstretch_rejections = 0;
        double lacc, uacc;
        while (i < NMC) {
            skip_trial = false;
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    if (!p0.linked()) {
                        overstretch_rejections = overstretch_rejections + 1.0;
                    }
                }
                else {
                    lacc = p0.link_acceptance(1);
                    if (lacc != 0.0) {
                        link_acceptances.push_back(lacc);
                    }
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
                uacc = p0.link_acceptance(0);
                if (skip_trial == false && uacc != 0.0) {
                    unlink_acceptances.push_back(p0.link_acceptance(0));
                }
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {

                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    linker_weights_l.push_back(p0.get_hairpin_weight(true));
                    helix_weights_l.push_back(p0.get_helix_weight());

                    old_linker_weights.push_back(p0.get_hairpin_weight(false));
                    new_linker_weights.push_back(p0.get_hairpin_weight(true));

                    if (l != 0) {
                        old_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, false));
                        new_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, true));
                    }


                }
            }

            i++;
        }

        linker_weights[l] = sum_of_elements(linker_weights_l) / static_cast<double>(linker_weights_l.size());
        helix_weights[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[l] = overstretch_rejections;


        double_vector_to_txt(path + "link_acc" + std::to_string(l) + ".txt", link_acceptances);
        double_vector_to_txt(path + "unlink_acc" + std::to_string(l) + ".txt", unlink_acceptances);


        double_vector_to_txt(path + "old_linker_weights_dangle=" + std::to_string(l) + ".txt", old_linker_weights);
        double_vector_to_txt(path + "new_linker_weights_dangle=" + std::to_string(l) + ".txt", new_linker_weights);

        double_vector_to_txt(path + "old_dangle_weights_dangle=" + std::to_string(l) + ".txt", old_dangle_weights);
        double_vector_to_txt(path + "new_dangle_weights_dangle=" + std::to_string(l) + ".txt", new_dangle_weights);

        double_vector_to_txt(path + "old_weights_dangle=" + std::to_string(l) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weights_dangle=" + std::to_string(l) + ".txt", new_weights);

    }
    double_vector_to_txt(path + "prefactors.txt", prefactors);
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);
    double_vector_to_txt(path + "linker_weights.txt", linker_weights);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);

}

void varying_helix_length(int N_monomers, int NMC, std::string path)
{
    int i{ 0 },l{0};
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    int helix_max_length{ 20 };

    int N_0{N_monomers};


    bool skip_trial{ false };
    std::vector<double> old_weights, new_weights;
    std::vector<double> linker_weights(helix_max_length-3), helix_weights(helix_max_length - 3), 
        linker_weights_l, helix_weights_l, physical_rejections(helix_max_length - 3),
            old_linker_weights, new_linker_weights, link_probability(helix_max_length-3);


    std::vector<double> link_probabilities(helix_max_length - 3), prefactors(helix_max_length-3);
    double link_prob;
    //p0->set_link_acceptance_prefactor(1);

    std::vector<int> helix;
    int N;

    while (l < helix_max_length-3) {
        if (l == 14-3) {
            std::cout << "stop" << std::endl;
        }


        std::cout << "considering helices of length" << 3 + l << std::endl;
        helix = { 0,2 + l,N_0 + 2 * l - 3 - l,N_0 + 2 * l - 1 };
        N = N_0 + 2 * l;

        polymer p0(N, true);
        p0.set_search_results(helix);

        prefactors[l] = get_prefactor(helix, N, NMC, 100);
        prefactors[l] = correlation_length(prefactors[l], N, helix, NMC, std::to_string(l));
        p0.set_link_acceptance_prefactor(prefactors[l]);

        i = 0;
        bound_counts = 0;
        unbound_counts = 0;


        double overstretch_rejections{ 0 };
        std::vector<double> link_acceptances, unlink_acceptances;


        helix_weights_l.clear();
        old_linker_weights.clear();
        new_linker_weights.clear();
        old_weights.clear();
        new_weights.clear();
        link_acceptances.clear();
        unlink_acceptances.clear();
        overstretch_rejections = 0;

        double lacc, uacc;
        while (i < NMC) {

            skip_trial = false;
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    if (!p0.linked()) {
                        overstretch_rejections = overstretch_rejections + 1.0;
                    }
                }
                else {
                    lacc = p0.link_acceptance(1);
                    if (lacc != 0.0) {
                        link_acceptances.push_back(lacc);
                    }
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
                uacc = p0.link_acceptance(0);
                if (skip_trial == false && uacc != 0.0) {
                    unlink_acceptances.push_back(p0.link_acceptance(0));
                }
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {

                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    helix_weights_l.push_back(p0.get_helix_weight());

                    old_linker_weights.push_back(p0.get_hairpin_weight(false));
                    new_linker_weights.push_back(p0.get_hairpin_weight(true));

                }
            }

            i++;
        }
        helix_weights[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[l] = overstretch_rejections;

        double_vector_to_txt(path + "link_acc_extra_bp=" + std::to_string(l) + ".txt", link_acceptances);
        double_vector_to_txt(path + "unlink_acc_extra_bp=" + std::to_string(l) + ".txt", unlink_acceptances);
        double_vector_to_txt(path + "old_linker_weights_extra_bp=" + std::to_string(l) + ".txt", old_linker_weights);
        double_vector_to_txt(path + "new_linker_weights_extra_bp=" + std::to_string(l) + ".txt", new_linker_weights);
        double_vector_to_txt(path + "old_weights_extra_bp=" + std::to_string(l) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weights_extra_bp=" + std::to_string(l) + ".txt", new_weights);

        l++;


    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

    double_vector_to_txt(path + "prefactors.txt", prefactors);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);



    print_1d_doub_vec(link_probabilities);
}

//void vary_dangle_sim(polymer* p0, int NMC, std::string filename)
//{
//
//    int max_dangle_length{ 20 };
//
//    //data 
//    std::vector<double> helix_weights(max_dangle_length), helix_weights_l(NMC), link_probability(max_dangle_length), link_acceptances, unlink_acceptances;
//
//    double overstretch_rejections{0};
//
//    bool skip_trial{ false };
//
//    int i{ 0 };
//    double u_branch;
//    int bound_counts{ 0 }, unbound_counts{ 0 };
//    int N{ p0->get_num_monomers() };
//    std::vector<double> old_weights, new_weights;
//    p0->force_unbound_state();
//
//
//
//    for (int l{ 0 }; l < max_dangle_length; l++) {
//        std::cout << "dangle length " << l << std::endl;
//        //p0->set_search_results({ 0 + l,2 + l,N - 3,N - 1 });
//        p0->set_search_results({ 0,2,N - 3-l,N - 1-l });
//
//        p0->force_unbound_state();//always begin in unlinked state
//        i = 0;
//        bound_counts = 0;
//        unbound_counts = 0;
//
//        old_weights.clear();
//        new_weights.clear();
//        link_acceptances.clear();
//        unlink_acceptances.clear();
//        overstretch_rejections = 0;
//
//        while (i < NMC) {
//            std::cout << i << "out of "<< NMC << std::endl;
//            u_branch = rand2(0, 1);
//            if (u_branch < 0.5) {// link_branch
//                std::cout << "link branch" << std::endl;
//                p0->link_move(skip_trial);
//                if (skip_trial == true) {
//                    overstretch_rejections = overstretch_rejections + 1.0;
//                }
//                else {
//                    link_acceptances.push_back(p0->link_acceptance(1));
//                }
//            }
//            else {
//                std::cout << "unlink branch" << std::endl;
//
//                p0->unlink_move(skip_trial);
//                if (skip_trial == false) {
//                    unlink_acceptances.push_back(p0->link_acceptance(0));
//                }
//            }
//
//            if (p0->linked() == true) {
//                bound_counts++;
//            }
//            else {
//                unbound_counts++;
//            }
//            if (skip_trial == false) {
//                double old_w{ p0->get_old_config_weight() }, new_w{ p0->get_new_config_weight() };
//
//                if (old_w != 0 && new_w != 0) {
//                    old_weights.push_back(p0->get_old_config_weight());
//                    new_weights.push_back(p0->get_new_config_weight());
//                }
//            }
//
//            helix_weights_l[i] = p0->get_helix_weight();
//            i++;
//        }
//        helix_weights[l] = sum_of_elements(helix_weights_l) / NMC;
//        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
//
//        double_vector_to_txt("Results/varying dangle sims/overstretch_rejections" + std::to_string(l) + ".txt", { overstretch_rejections });
//        double_vector_to_txt("Results/varying dangle sims/link_acc" + std::to_string(l) + ".txt", link_acceptances);
//        double_vector_to_txt("Results/varying dangle sims/unlink_acc" + std::to_string(l) + ".txt", unlink_acceptances);
//
//
//        double_vector_to_txt("Results/varying dangle sims/old_weights_dangle=" + std::to_string(l) + ".txt", old_weights);
//        double_vector_to_txt("Results/varying dangle sims/new_weights_dangle=" + std::to_string(l) + ".txt", new_weights);
//
//    }
//    print_1d_doub_vec(helix_weights);
//    double_vector_to_txt("Results/varying dangle sims/varying_dangle_helix_weights.txt", helix_weights);
//    double_vector_to_txt("Results/varying dangle sims/varying_dangle_link_probability.txt", link_probability);
//    //std::cout << "Unlink state: " << unbound_counts << std::endl;
//    //std::cout << "Link state: " << bound_counts << std::endl;
//
//}

//void vary_dangle_sim(double acceptance_prefactor, int N_monomers, int NMC, const std::vector<int>& search_result, std::string path)
//{
//    int max_dangle_length{ 20 };
//
//    //data 
//    std::vector<double> link_probability(max_dangle_length), link_acceptances, unlink_acceptances;
//
//    std::vector<double> linker_weights(max_dangle_length), helix_weights(max_dangle_length), linker_weights_l, helix_weights_l, physical_rejections(max_dangle_length);
//    double overstretch_rejections{ 0 };
//
//    bool skip_trial{ false };
//
//
//    int i{ 0 };
//    double u_branch;
//    int bound_counts{ 0 }, unbound_counts{ 0 };
//    int N{ N_monomers };
//    std::vector<double> old_weights, new_weights;
//    //p0->force_unbound_state();
//
//
//
//    for (int l{ 0 }; l < max_dangle_length; l++) {
//        polymer p0(N_monomers, true);
//        p0.set_link_acceptance_prefactor(acceptance_prefactor);
//        std::cout << "dangle length " << l << std::endl;
//        p0.set_search_results({ 0 + l,2 + l,N - 3,N - 1 });
//        //p0.set_search_results(search_result);
//
//        //p0->force_unbound_state();//always begin in unlinked state
//        i = 0;
//        bound_counts = 0;
//        unbound_counts = 0;
//
//        helix_weights_l.clear();
//        linker_weights_l.clear();
//        old_weights.clear();
//        new_weights.clear();
//        link_acceptances.clear();
//        unlink_acceptances.clear();
//        overstretch_rejections = 0;
//        double lacc, uacc;
//        while (i < NMC) {
//            skip_trial = false;
//            std::cout << i << "out of " << NMC << std::endl;
//            u_branch = rand2(0, 1);
//            if (u_branch < 0.5) {// link_branch
//                std::cout << "link branch" << std::endl;
//                p0.link_move(skip_trial);
//                if (skip_trial == true) {
//                    if (!p0.linked()) {
//                        overstretch_rejections = overstretch_rejections + 1.0;
//                    }
//                }
//                else {
//                    lacc = p0.link_acceptance(1);
//                    if (lacc != 0.0) {
//                        link_acceptances.push_back(lacc);
//                    }
//                }
//            }
//            else {
//                std::cout << "unlink branch" << std::endl;
//
//                p0.unlink_move(skip_trial);
//                uacc = p0.link_acceptance(0);
//                if (skip_trial == false && uacc != 0.0) {
//                    unlink_acceptances.push_back(p0.link_acceptance(0));
//                }
//            }
//
//            if (p0.linked() == true) {
//                bound_counts++;
//            }
//            else {
//                unbound_counts++;
//            }
//            if (skip_trial == false) {
//                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };
//
//                if (old_w != 0 && new_w != 0) {
//
//                    old_weights.push_back(p0.get_old_config_weight());
//                    new_weights.push_back(p0.get_new_config_weight());
//
//                    linker_weights_l.push_back(p0.get_hairpin_weight(true));
//                    helix_weights_l.push_back(p0.get_helix_weight());
//
//                }
//            }
//
//            i++;
//        }
//
//        linker_weights[l] = sum_of_elements(linker_weights_l) / static_cast<double>(linker_weights_l.size());
//        helix_weights[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
//        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
//        physical_rejections[l] = overstretch_rejections;
//
//
//        double_vector_to_txt(path + "link_acc" + std::to_string(l) + ".txt", link_acceptances);
//        double_vector_to_txt(path + "unlink_acc" + std::to_string(l) + ".txt", unlink_acceptances);
//
//
//        double_vector_to_txt(path + "old_weights_dangle=" + std::to_string(l) + ".txt", old_weights);
//        double_vector_to_txt(path + "new_weights_dangle=" + std::to_string(l) + ".txt", new_weights);
//
//    }
//
//    double_vector_to_txt(path + "helix_weights.txt", helix_weights);
//    double_vector_to_txt(path + "linker_weights.txt", linker_weights);
//    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
//    double_vector_to_txt(path + "link_probabilities.txt", link_probability);
//
//}

void vary_dangle_sim(int N_monomers, int NMC, std::string path)
{
    int max_dangle_length{ 20 };

    //data 
    std::vector<double> link_probability(max_dangle_length), link_acceptances, unlink_acceptances, prefactors(max_dangle_length);

    std::vector<double> linker_weights(max_dangle_length), helix_weights(max_dangle_length), 
        linker_weights_l, helix_weights_l, physical_rejections(max_dangle_length),
        old_linker_weights, new_linker_weights, old_dangle_weights, new_dangle_weights;
    double overstretch_rejections{ 0 };

    bool skip_trial{ false };


    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    std::vector<double> old_weights, new_weights;
    //p0->force_unbound_state();

    int N0{ N_monomers }, N{N0};
    std::vector<int> helix(4);

    for (int l{ 0 }; l < max_dangle_length; l++) {

        // the tail is going to get longer such that the total number of monomers in the polymer increases each time.
        N = N0 + 2*l;
        polymer p0(N, true);
        helix = { 0 + l, 2 + l, N - 3-l, N - 1-l };
        p0.set_search_results(helix);



        prefactors[l] = get_prefactor(helix, N, NMC, 100);
        prefactors[l] = correlation_length(prefactors[l], N, helix, NMC, std::to_string(l));
        p0.set_link_acceptance_prefactor(prefactors[l]);

        std::cout << "dangle length " << l << std::endl;

        i = 0;
        bound_counts = 0;
        unbound_counts = 0;

        helix_weights_l.clear();
        linker_weights_l.clear();
        old_weights.clear();
        new_weights.clear();
        link_acceptances.clear();
        unlink_acceptances.clear();
        overstretch_rejections = 0;
        double lacc, uacc;
        while (i < NMC) {
            skip_trial = false;
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    if (!p0.linked()) {
                        overstretch_rejections = overstretch_rejections + 1.0;
                    }
                }
                else {
                    lacc = p0.link_acceptance(1);
                    if (lacc != 0.0) {
                        link_acceptances.push_back(lacc);
                    }
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
                uacc = p0.link_acceptance(0);
                if (skip_trial == false && uacc != 0.0) {
                    unlink_acceptances.push_back(p0.link_acceptance(0));
                }
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {

                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    linker_weights_l.push_back(p0.get_hairpin_weight(true));
                    helix_weights_l.push_back(p0.get_helix_weight());

                    old_linker_weights.push_back(p0.get_hairpin_weight(false));
                    new_linker_weights.push_back(p0.get_hairpin_weight(true));

                    if (l != 0) {
                        old_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, false));
                        new_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, true));
                    }


                }
            }

            i++;
        }

        linker_weights[l] = sum_of_elements(linker_weights_l) / static_cast<double>(linker_weights_l.size());
        helix_weights[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[l] = overstretch_rejections;


        double_vector_to_txt(path + "link_acc" + std::to_string(l) + ".txt", link_acceptances);
        double_vector_to_txt(path + "unlink_acc" + std::to_string(l) + ".txt", unlink_acceptances);


        double_vector_to_txt(path + "old_linker_weights_dangle=" + std::to_string(l) + ".txt", old_linker_weights);
        double_vector_to_txt(path + "new_linker_weights_dangle=" + std::to_string(l) + ".txt", new_linker_weights);

        double_vector_to_txt(path + "old_dangle_weights_dangle=" + std::to_string(l) + ".txt", old_dangle_weights);
        double_vector_to_txt(path + "new_dangle_weights_dangle=" + std::to_string(l) + ".txt", new_dangle_weights);

        double_vector_to_txt(path + "old_weights_dangle=" + std::to_string(l) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weights_dangle=" + std::to_string(l) + ".txt", new_weights);

    }
    double_vector_to_txt(path + "prefactors.txt", prefactors);
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);
    double_vector_to_txt(path + "linker_weights.txt", linker_weights);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);

}


void kissing_hairpin_swivel(polymer* p0, int NMC, bool swivel_allowed, std::string filename)
{
    double swivel_switch;
    swivel_allowed == true ? swivel_switch = 1 : swivel_switch = 2.0 / 3.0;
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    int kissing_hairpin_count{ 0 }, two_link_state{0}, one_link_state{0}, free_state{0};

    p0->force_unbound_state();

    while (i < NMC) {
        std::cout << i << std::endl;
        u_branch = rand2(0, swivel_switch);
        if (u_branch < 1.0/3.0) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0->link_move(m);
        }
        else if(u_branch <2.0/3.0) {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0->unlink_move(m);
        }
        else if (u_branch < 1) {
            p0->swivel_move();
        }

        if (p0->get_num_helices() == 0) {
            free_state++;
        }
        else if (p0->get_num_helices() == 1) {
            one_link_state++;
        }
        else if (p0->get_num_helices() == 2) {
            two_link_state++;
        }
        else if (p0->get_num_helices() == 3) {
            kissing_hairpin_count++;
        }
        i++;
    }
    std::cout << "out of main loop" << std::endl; 

    double P_kissing_hairpin{ static_cast<double>(kissing_hairpin_count) / static_cast<double>(NMC) }, 
    P_two_links{ static_cast<double>(two_link_state) / static_cast<double>(NMC) },
        P_one_link{ static_cast<double>(one_link_state) / static_cast<double>(NMC) },
        P_free{ static_cast<double>(free_state) / static_cast<double>(NMC) };
    std::string file_name;
    if (swivel_allowed == true) {
        file_name = "P_kissing_hairpin_ON.txt";
    }
    else {
        file_name = "P_kissing_hairpin_OFF.txt";
    }
    double_vector_to_txt("Results/" + file_name, {P_free, P_one_link, P_two_links,P_kissing_hairpin});

}

void old_new_weights_sim(polymer* p0, int NMC, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };


    std::vector<double> old_config_weights, new_config_weights;
    bool skip_trial{false};

    p0->force_unbound_state();

    while (i < NMC) {
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            p0->link_move(skip_trial);
        }
        else {

            p0->unlink_move(skip_trial);
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }
        if (skip_trial == false) {
            double old_w{ p0->get_old_config_weight() }, new_w{ p0->get_new_config_weight() };
            
            if (old_w != 0 && new_w != 0) {
                old_config_weights.push_back(p0->get_old_config_weight());
                new_config_weights.push_back(p0->get_new_config_weight());
            }
        }
        i++;
    }

    double_vector_to_txt("Results/old_vs_new_old_weights.txt", old_config_weights);
    double_vector_to_txt("Results/old_vs_new_new_weights.txt", new_config_weights);



}

void vary_WCA_potential(double prefactor, int N_monomers, int NMC, const std::vector<int>& search_result, std::string path)
{
    int i{ 0 }, k{ 0 };
    double n_runs{ 10.0 };
    std::vector<double> epsilons{ generate_linear_array(0,2,10) };

    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    bool skip_trial{ false };
    std::vector<double> physical_rejections(n_runs);
    int overstretch_rejections{ 0 };

    std::vector<double> old_weights, new_weights, helix_weights_l, linker_weights_l;
    std::vector<double> helix_weights(n_runs), linker_weights(n_runs), link_probability(n_runs), avg_linker_weight(n_runs);

    double current_epsilon;
    while (k < n_runs) {
        polymer p0(N_monomers, true);
        p0.set_link_acceptance_prefactor(prefactor);
        p0.set_search_results(search_result);

        current_epsilon = epsilons[k];
        set_WCA_parameter(current_epsilon);


        // reset temp data for new round
        i = 0;
        bound_counts = 0;
        unbound_counts = 0;
        linker_weights_l.clear();
        helix_weights_l.clear();
        old_weights.clear();
        new_weights.clear();
        overstretch_rejections = 0;

        while (i < NMC) {
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    overstretch_rejections = overstretch_rejections + 1.0;
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {
                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    linker_weights_l.push_back(p0.get_hairpin_weight(true));
                    helix_weights_l.push_back(p0.get_helix_weight());
                }
            }

            i++;
        }
        linker_weights[k] = sum_of_elements(linker_weights_l) / static_cast<double>(linker_weights_l.size());
        helix_weights[k] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[k] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[k] = overstretch_rejections;

        double_vector_to_txt(path + "old_weightsEPS=" + std::to_string(current_epsilon) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weightsEPS=" + std::to_string(current_epsilon) + ".txt", new_weights);


        k++;
    }
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);
    double_vector_to_txt(path + "linker_weights.txt", linker_weights);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);

}

void testing_convergence(int NMC, std::string path)
{

    std::vector<int> s;

    int max_linker_length{ 15 }, N0{8},N{N0};
    int update_rate{ 100 };
    double delta{ helix_separation() };

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };


    std::vector<double> old_config_weights, new_config_weights, linker_weights, helix_weights;
    bool skip_trial{ false };

    std::vector<double> ideal_pdf(max_linker_length), prefactors(max_linker_length), p_linked(max_linker_length);
    

    for (int l{ 0 }; l < max_linker_length; l++) {
        N = N0 + l;
        
        polymer p0(N, true);
        s = { 0,2,5 + l,7 + l };


        ideal_pdf[l] = ideal_chain_pdf(delta, l + 3);
        prefactors[l] = get_prefactor(s, N, NMC, 200);
        std::cout << " out of prefactor " << std::endl;
        p0.set_link_acceptance_prefactor(prefactors[l]);
        p0.set_search_results(s);

        i = 0;
        bound_counts = 0;
        unbound_counts = 0;

        old_config_weights.clear();
        new_config_weights.clear();
        helix_weights.clear();
        linker_weights.clear();
        while (i < NMC) {
            skip_trial = false;
            std::cout << i << "out of " << NMC << std::endl;

            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                p0.link_move(skip_trial);
            }
            else {

                p0.unlink_move(skip_trial);
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };
                if (old_w != 1 || new_w != 1) {
                    p0.get_old_config_weight();
                    p0.get_new_config_weight();
                }
                if (old_w != 0 && new_w != 0) {
                    old_config_weights.push_back(p0.get_old_config_weight());
                    new_config_weights.push_back(p0.get_new_config_weight());

                    linker_weights.push_back(p0.get_hairpin_weight(true));
                    helix_weights.push_back(p0.get_helix_weight());

                }
            }
            i++;





        }
        p_linked[l] = static_cast<double>(bound_counts) / (static_cast<double>(bound_counts) + static_cast<double>(unbound_counts));


        double_vector_to_txt(path + "helix_weights" + std::to_string(l) + ".txt", helix_weights);
        double_vector_to_txt(path + "linker_weights" + std::to_string(l) + ".txt", linker_weights);
        double_vector_to_txt(path + "old_weights" + std::to_string(l) + ".txt", old_config_weights);
        double_vector_to_txt(path + "new_weights" + std::to_string(l) + ".txt", new_config_weights);

    }
    double_vector_to_txt(path + "p_link.txt", p_linked);
    double_vector_to_txt(path + "ideal_chain_pdf.txt", ideal_pdf);
    double_vector_to_txt(path + "prefactors.txt", prefactors);


}

