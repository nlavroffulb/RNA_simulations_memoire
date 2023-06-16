#include "polymer_class.h"
#include <algorithm>
#include <ctime>


//polymer constructors
//**************************************************************************************//
//**************************************************************************************//
//**************************************************************************************//
polymer::polymer(int num) {
    N = num;
    std::vector<std::vector<double>> positions(N);
    std::vector<double> weights(N), energies(N);
    std::vector<double> start{ 0.0,0.0,0.0 };
    positions = rosenbluth_random_walk(start, N - 1, excluded_volume, energies, weights);
    std::vector<monomer*> M(N);
    for (int i{ 0 }; i < N; i++) {
        M[i] = new monomer(positions[i], i, random_base());
        M[i]->left_right_link[1] = positions.size() - 1;
    }
    chain = M;
    search_results = structure_search();

    std::vector<int> limit{ 0,N - 1 };
    old_weights = new rosenbluth_growth(weights, energies);
    new_weights = new rosenbluth_growth(weights, energies);

    new_config_positions = positions;
}

polymer::polymer(int num, bool rosenbluth, int struct_size)
{
    N = num;
    std::vector<std::vector<double>> positions(N);
    std::vector<double> weights(N), energies(N);
    std::vector<double> start{ 0.0,0.0,0.0 };
    excluded_volume.push_back(start);
    if (rosenbluth == true) {
        positions = rosenbluth_random_walk(start, N - 1, excluded_volume, energies, weights);
        energies.insert(energies.begin(), 0);
        weights.insert(weights.begin(), 1);
    }
    
    else if (rosenbluth == false) {
        positions = random_walk(start, N - 1, excluded_volume);
    }
    std::vector<monomer*> M(N);
    for (int i{ 0 }; i < N; i++) {
        M[i] = new monomer(positions[i], i, random_base());
        M[i]->left_right_link[1] = positions.size() - 1;
    }
    chain = M;
    
    search_results = structure_search(struct_size);

    rosenbluth_switch = rosenbluth;
    std::vector<int> limit{ 0,N - 1 };
    if (rosenbluth_switch == true) {
        old_weights = new rosenbluth_growth(weights, energies);
        new_weights = new rosenbluth_growth(weights, energies);

    }
    else {
        std::fill(weights.begin(), weights.end(), 1.0);
        std::fill(energies.begin(), energies.end(), 0.0);

        old_weights = new rosenbluth_growth(weights, energies);
        new_weights = new rosenbluth_growth(weights, energies);


    }

    new_config_positions = positions;
    old_config_positions = positions;
}

polymer::polymer(int num, bool rosenbluth, std::vector<std::string>& sequence, int struct_size)
{
    N = num;
    std::vector<std::vector<double>> positions(N);
    std::vector<double> weights(N), energies(N);
    std::vector<double> start{ 0.0,0.0,0.0 };
    excluded_volume.push_back(start);
    if (rosenbluth == true) {
        positions = rosenbluth_random_walk(start, N - 1, excluded_volume, energies, weights);
        energies.insert(energies.begin(), 0);
        weights.insert(weights.begin(), 1);
    }

    else if (rosenbluth == false) {
        positions = random_walk(start, N - 1, excluded_volume);
    }
    std::vector<monomer*> M(N);
    for (int i{ 0 }; i < N; i++) {
        M[i] = new monomer(positions[i], i, sequence[i]);
        M[i]->left_right_link[1] = positions.size() - 1;
    }
    chain = M;

    search_results = structure_search(struct_size);

    rosenbluth_switch = rosenbluth;
    std::vector<int> limit{ 0,N - 1 };
    if (rosenbluth_switch == true) {
        old_weights = new rosenbluth_growth(weights, energies);
        new_weights = new rosenbluth_growth(weights, energies);

    }
    else {
        std::fill(weights.begin(), weights.end(), 1.0);
        std::fill(energies.begin(), energies.end(), 0.0);

        old_weights = new rosenbluth_growth(weights, energies);
        new_weights = new rosenbluth_growth(weights, energies);


    }

    new_config_positions = positions;
    old_config_positions = positions;

}



polymer::~polymer()
{
    for (auto i : chain) {
        delete i;
        //std::cout << chain.size() << std::endl;
    }
    std::cout << "Polymer destructor called." << std::endl;
}


//**************************************************************************************//
//**************************************************************************************//
//**************************************************************************************//


//general functions of the polymer and given monomer in the chain



std::vector<std::vector<double>> polymer::generate_vector_of_monomer_positions()
{
    std::vector<std::vector<double>> temp;
    for (int i{ 0 }; i < chain.size(); i++) {
        temp.push_back(chain[i]->get_position());
    }

    return temp;
}
void polymer::print_monomer_positions()
{
    for (int i{ 0 }; i < chain.size(); i++) {
        std::cout << i << " ";

        std::cout << "{ ";

        for (int j{ 0 }; j < 3; j++) {
            std::cout << generate_vector_of_monomer_positions()[i][j] << " ";
        }
        std::cout << "}" << std::endl;
    }
}

void polymer::print_bases()
{
    for (auto m : chain) {
        std::cout << m->get_base() << " ";
    }/*std::cout << std::endl;*/

    std::cout << "\n";
}
void polymer::print_bases1()
{
    for (auto m : chain) {
        std::cout << m->get_id()<< " "<< m->get_base() << std::endl;
    }/*std::cout << std::endl;*/

    std::cout << "\n";
}


void polymer::output_for_ovito4(std::string filename) {

    std::vector<int> linked_monomers;
    get_linked_monomers(linked_monomers);
    //std::vector<int> s_temp(4);
    //for (auto s : helix_list) {
    //    s_temp = s->get_monomers();
    //    unpack_region(s_temp);
    //    for (auto i : s_temp) {
    //        linked_monomers.push_back(i);
    //    }
    //}
    std::string file_name{ "polymer " + std::to_string(rand2(0,1)) + ".txt" };

    if (filename != "") {
        file_name = filename + ".txt";
    }

    std::ofstream ovito_output(file_name);
    std::vector<std::vector<double>> polymer{ generate_vector_of_monomer_positions() };
    ovito_output << polymer.size() << "\n\n";
    for (int i{ 0 }; i < polymer.size(); i++) {
        int c{ 0 }, d{ 0 };

        std::string line{ "" };
        for (int j{ 0 }; j < linked_monomers.size(); j++) {
            if (i == linked_monomers[j]) {
                c++;

            }
        }
        //if (i > linked_monomers[linked_monomers.size() / 2 - 1] && i < linked_monomers[structure.size() / 2]) {
        //    //std::cout << "working" << std::endl;
        //    line += "Mg";
        //}
        if (i == 0) {
            line += "As";
            //continue;
        }
        //else if (i == chain.size() - 2) {
        //    line += "Mg";
        //}
        else if (i == chain.size() - 1) {
            line += "As";
        }
        else if (c > 0) {
            line += "He";
        }
        else if (i % 2 == 0) {
            line += "C";
        }
        else if (i % 2 != 0) {
            line += "O";
        }
        for (int j{ 0 }; j < polymer[i].size(); j++) {
            line += " " + std::to_string(polymer[i][j]);
        }
        if (i == polymer.size() - 1) {
            //std::cout << "this line ran";
            ovito_output << line << std::endl;;
        }
        else {
            ovito_output << line << "\n";
        }
    }

    ovito_output.close();

}
void polymer::ovito_bonds(std::string filename)
{
    std::vector<int> linked_monomers;
    get_linked_monomers(linked_monomers);
    std::string file_name{ "polymer " + std::to_string(rand2(0,1)) + ".txt" };

    if (filename != "") {
        file_name = filename + ".txt";
    }

    // this block adds all the particles that we want
    std::ofstream ovito_output(file_name);
    std::vector<std::vector<double>> polymer{ generate_vector_of_monomer_positions() };
    int np{ 20 };// the number of particles we'll add to represent a bond

    int n_particles{0};

    std::vector<int> dh(4);
    for (auto i : helix_list) {// count the number of bond particles in our double helix structures
        dh = i->get_monomers();
        n_particles += np * (abs(dh[1] - dh[0]) + 1);
    }
    n_particles += polymer.size(); // the actual particles
    n_particles += np * (polymer.size() - 1); // the bond particles due to bonds between consecutive monomers
    ovito_output << n_particles<< "\n\n";


    for (int i{ 0 }; i < polymer.size(); i++) {
        int c{ 0 }, d{ 0 };

        std::string line{ "" };
        for (int j{ 0 }; j < linked_monomers.size(); j++) {
            if (i == linked_monomers[j]) {
                c++;

            }
        }
        if (i == 0) {
            line += "As";
        }
        else if (i == chain.size() - 1) {
            line += "As";
        }
        else if (c > 0) {
            line += "He";
        }
        else {
            line += "O";
        }
        for (int j{ 0 }; j < polymer[i].size(); j++) {
            line += " " + std::to_string(polymer[i][j]);
        }
        ovito_output << line << "\n";

    }

    // this block will add particles to show bonds. essentially 2 small particles in between each consecutive monomer as well as within structures.

    //consecutive monomers
    for (int i = 0; i < polymer.size(); i++) {

        if (i != polymer.size() - 1) {

            int j{ i + 1 };

            // calculate the vector which goes from one to the other. 
            // then divide the line into np + 1.
            std::vector<double> v{ vector_subtraction(polymer[j],polymer[i]) };
            double modv{ vector_modulus(v) };
            modv = modv / (np + 1);
            std::vector<double> temp(3), bond_pos(3);
            for (int k{ 1 }; k <= np; k++) {
                std::string line{ "" };
                line += "Mg";

                temp = multiplication_by_scalar(k*modv, v);
                bond_pos = vector_addition(polymer[i], temp);

                for (int l{ 0 }; l < bond_pos.size(); l++) {
                    line += " " + std::to_string(bond_pos[l]);
                }

                //if(i == polymer.size()-2 && k == np){
                //    ovito_output << line << std::endl;

                //}
                //else {
                //    ovito_output << line << "\n";
                //    
                //}
                ovito_output << line << "\n";


            }
        }
        //else {
        //    std::string line{ "" };
        //    ovito_output << line << "\n";
        //}
    }
    //bonds in the helix
    //"Ti"
    for (int i{ 0 }; i < helix_list.size();i++) {
        std::vector<int> monomers{ helix_list[i]->get_monomers() };
        unpack_region(monomers);

        // bonds between pairs
        for (int k{ 0 }; k < monomers.size() /2 ; k++){
            int m1{ monomers[0] + k }, m2{monomers.back() - k};

            // calculate the vector which goes from one to the other. 
            // then divide the line into np + 1.
            std::vector<double> v{ vector_subtraction(polymer[m2],polymer[m1]) };
            double modv{ vector_modulus(v) };
            modv = modv / (np + 1);
            std::vector<double> temp(3), bond_pos(3);

            for (int l{ 1 }; l <= np; l++) {
                std::string line{ "" };
                line += "Ti";

                temp = multiplication_by_scalar(l * modv, v);
                bond_pos = vector_addition(polymer[m1], temp);

                for (int p{ 0 }; p < bond_pos.size(); p++) {
                    line += " " + std::to_string(bond_pos[p]);
                }

                if (i == helix_list.size() - 1 && k == -1 + monomers.size() / 2 && l == np) {
                    ovito_output << line << std::endl;
                }
                else {
                    ovito_output << line << "\n";
                }
            }
        }
    }
    ovito_output.close();


}
std::vector<std::string> polymer::get_sequence()
{
    std::vector<std::string> sequence(N);
    for (int i = 0; i < chain.size(); i++)
    {
        sequence[i] = chain[i]->get_base();
    }
    return sequence;
}
std::vector<std::vector<double>> polymer::get_subpositions(std::vector<int> limit)
{
    int kappa;
    limit[1] < limit[0] ? kappa = -1 : kappa = 1;
    std::vector<std::vector<double>> temp(abs(limit[1] - limit[0]) + 1);
    for (int i{ 0 }; i < temp.size(); i++) {
        temp[i] = chain[limit[0]+kappa*i]->get_position();
    }
    return temp;
}
std::vector<std::string> particle_types{ "As","Gd","Hf","Mg","K", "P"};

monomer* polymer::operator[](int i) {
    return chain[i];
}
//**************************************************************************************//
//**************************************************************************************//
//**************************************************************************************//


//**************************************************************************************//
//**************************************************************************************//
//**************************************************************************************//

void polymer::get_linked_monomers(std::vector<int>& links) {
    std::vector<int> s(4);
    int limits;

    if (helix_list.size() > 0) {
        std::vector<int> linked_monomers;
        for (auto i : helix_list) {
            s = i->get_monomers();
            unpack_region(s);

            for (auto j : s) {
                linked_monomers.push_back(j);
            }
            //limits = s[1] - s[0];
            //for (int j{ 0 }; j < limits; j++) {
            //    linked_monomers.push_back(s[0] + j);
            //    linked_monomers.push_back(s[2] + j);
            //}
        }
        links = linked_monomers;
    }
}
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//


std::vector<std::vector<int>> polymer::structure_search(int n)//default argument is n=3. it also returns the search results
// always as a 4 element int vector delimiting the compatible region 1 [a,b] and 2 [c,d] => search result = [a,b,c,d]
{

    std::vector<std::string> bases(chain.size());
    for (int b{ 0 }; b < chain.size(); b++) {
        bases[b] = chain[b]->get_base();
        //bases.push_back(chain[b]->get_base());
    }

    std::vector<std::string> regions_r;
    //for (int i{ 0 }; i < r_max-1; i++) {
    for (int i{ 0 }; i < bases.size() - 1; i++) {
        if (chain.size() - i <= n - 1) {//when we get towards the end of the chain, e.g r=2 for U A G C, once we get to C we cannot form a group of 2 so break.
            break;
        }
        else {
            std::string region{ "" };
            for (int j{ 0 }; j < n; j++) {

                region += bases[i + j];

            }
            regions_r.push_back(region);

        }
    }
    std::vector<std::vector<int>> search_results_r;//we store the compatible regions for a given group size r.

    for (int k{ 0 }; k < regions_r.size() - 1; k++) {//loop through regions. k is the region we compare all the others to(l).
        for (int l{ n }; l < regions_r.size(); l++) {

            //we don't need to compare compatibility for regions l which share bases with region k.
            //so we can skip that part of the loop with the continue statement.
            if (abs(l - k) < n || l < k) {
                continue;
            }

            //here we save the corresponding position (integer order in the chain, not euclidian position) of compatible bases
            if (regions_compatible(regions_r[k], regions_r[l])) {
                std::vector<int> format_result(4);

                format_result[0] = k;
                format_result[1] = k + n-1;
                format_result[2] = l;
                format_result[3] = l + n -1;

                search_results_r.push_back(format_result);
            }
        }
    }
    //if (search_results_r.size() > 0) {
    //    std::cout << "Results for region of length " << n << std::endl;
    //    for (auto v : search_results_r) {
    //        for (auto e : v) {
    //            std::cout << e << " ";
    //        }std::cout << std::endl;
    //    }//std::cout << std::endl;

    //}
    //print_2d_int_vec(search_results_r);
    //std::cout << "regions found " << search_results_r.size() << std::endl;
    return search_results_r;
}

void polymer::set_search_results(int n)
{
}

void polymer::sample_unzip_region(std::vector<int>& region_to_unzip, int side) {
    std::vector<int> contraction{ region_to_unzip };

    if (side == 0) {
        contraction[0] = contraction[0] + 1;
        contraction[3] = contraction[3] - 1;
    }
    else if (side == 1) {
        contraction[1] = contraction[1] - 1;
        contraction[2] = contraction[2] + 1;
    }


    region_to_unzip = contraction;


}
void polymer::sample_zip_region(std::vector<int>& region_to_zip, std::vector<int>& xtend_limits, int side) {
    std::vector<int> extension{ region_to_zip };

    if (side == 0) {
        extension[0] = extension[0] - 1;
        extension[3] = extension[3] + 1;
    }
    else if (side == 1) {
        extension[1] = extension[1] + 1;
        extension[2] = extension[2] - 1;
    }
    if (extension[0] < xtend_limits[0] || extension[1] > xtend_limits[1]) {
        region_to_zip = {};
    }
    else {
        region_to_zip = extension;
    }

    

}
void polymer::sample_link_region(std::vector<int>& link, int& alpha, int& beta, int &ss_index)
{
    //std::cout << "there are " << helix_list.size() << " linked structures" << std::endl;

    if (search_results.size() == 0) {
        link = {};
        return;
    }

    //// we're only considering regions of length 3 for the moment. the default argument of search results function is for r=3.
    int region_index0{ static_cast<int>(rand2(0,search_results.size())) };
    ss_index = region_index0;
    link = search_results[region_index0];

    // we only sample alpha and beta if there is more than 1 double helix structure already.
    if (helix_list.size() > 1) {
        alpha = static_cast<int>(rand2(0, 2));
        beta = static_cast<int>(rand2(0, 2));

    }
    else {
        alpha = 0;
        beta = 0;
    }
    std::cout << "alpha = " << alpha << " beta = " << beta << std::endl;

}

//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//

bool polymer::overlapping_monomers(std::vector<int> &reaction_region) {
    std::vector<int> linked_monomers;
    get_linked_monomers(linked_monomers);

    for (auto i : reaction_region) {
        if (std::any_of(linked_monomers.begin(), linked_monomers.end(), [i](int y) { return i == y; })) {
            //std::cout << "overlap of proposed structure with existing linked structure" << std::endl;
            return true;
        }
    }
    return false;

}


void polymer::sample_unlink_region(int &helix_index)
{
    std::cout << "there are " << helix_list.size() << " linked structures" << std::endl;
    std::vector<int> lm;
    get_linked_monomers(lm);
    print_1d_int_vec(lm);
    int s_i{ static_cast<int>(rand2(0,helix_list.size())) };
    if (helix_list[s_i]->extended()) {
        helix_index = -1;
    }else{
        helix_index = s_i;
    }
}

void polymer::unlink(int struct_index) {
    if (helix_list.size() == 0) {
        std::cout << "Cannot perform an unlink move" << std::endl;
    }
    else {
        //neighbouring_linkers();

        //pick structure to unlink and an origin within the structure to start the 'degrowth'
        //int s_i{ static_cast<int>(rand2(0,helix_list.size())) };
        std::vector<int> destructure{ helix_list[struct_index]->get_monomers() };

        // for this simulation, we don't need to sample alpha and beta
        int alpha{0}, beta{0};
        if (helix_list.size() !=1) {
            alpha = static_cast<int>(rand2(0, 2));
            beta = static_cast<int>(rand2(0, 2));

        }
        //std::cout << "alpha = " << alpha << " beta = " << beta << std::endl;

        std::vector<int> reaction_region{ destructure };
        int kappa{ static_cast<int>(pow(-1.0,alpha)) };
        int origin{ reaction_region[alpha * 3 + kappa * beta] }, origin_y{ reaction_region[(1 - alpha) * 3 - kappa * (1 - beta)] };

        // define a,b,c,d as the extremes of the structure to be destroyed. either a,b,c,d could be the origin of the unlink.
        int a, b, c, d;
        if (alpha == 1) {
            std::sort(reaction_region.begin(), reaction_region.end(), std::greater<>());

        }

        a = reaction_region[0];
        b = reaction_region[1];
        c = reaction_region[2];
        d = reaction_region[3];

        if (alpha == 1) {
            std::sort(reaction_region.begin(), reaction_region.end());

        }

        //define w,x,y,z as the leftmost linker, middle linkers and rightmost linkers respectively in the most
        //complicated case where there are structures on either side of the structure to be unlinked as well as in
        //the middle.
        int w, x, y, z;//sometimes some of these will be the same.
        w = chain[a]->left_right_link[alpha];
        x = chain[b]->left_right_link[1 - alpha];
        y = chain[c]->left_right_link[alpha];
        z = chain[d]->left_right_link[1 - alpha];

        //std::cout << "abcd " << a << " " << b << " " << c << " " << d << std::endl;
        //std::cout << "wxyz " << w << " " << x << " " << y << " " << z << std::endl;

        int region_size{ abs(b - a) + 1 };
        bool middle_struct{ true };
        if (abs(b - c) < abs(b - x)) {
            //x = c;
            y = b;
            middle_struct = false;
        }

        std::vector<std::vector<int>> degrowth_limits;

        if (beta == 0) {
            degrowth_limits.push_back({ a,x });
            if (middle_struct) {
                degrowth_limits.push_back({ y,z });
            }
        }
        else if (beta == 1) {
            if (middle_struct) {
                degrowth_limits.push_back({ w,b });
                degrowth_limits.push_back({ y,z });
            }
            else {
                degrowth_limits.push_back({ w,z });
            }

        }

        //std::cout << "degrowth limits" << std::endl;
        ////std::vector<std::vector<int>> temp;
        for (auto i : degrowth_limits) {
            std::cout << i[0] << " " << i[1] << std::endl;
        }
       
        update_excluded_volume(degrowth_limits);

        std::vector<std::vector<double>> s_x, s_y;
        std::vector<double> u, v;
        std::vector<double> o_h{ chain[origin]->get_position() };
        std::vector<std::vector<double>> rcentres(2);


        // grow the new structure (unlinked)
        if (rosenbluth_switch == false) {
            // need to get the old excluded volume interaction term for the helix.
            std::vector<std::vector<double>> helix_positions(2 * (abs(destructure[1] - destructure[0]) + 1));
            unpack_region(destructure);
            for (int i{ 0 }; i < destructure.size(); i++) {
                helix_positions[i] = chain[i]->get_position();
            }
            pack_region(destructure);
            helix_excluded_volume_interaction(helix_positions);

        }
        grow_limits(degrowth_limits, alpha, true);

        // grow the old structure (linked)
        std::vector<std::vector<int>> link_limits;
        link_growth_limits(destructure, alpha, beta, link_limits);
        update_excluded_volume(link_limits, destructure);
        // if we've turned on rosenbluth sampling then we rosenbluth sample the helix. the rosenbluth sample helix
        // returns the selected helix vectors by reference.
        if (rosenbluth_switch == true) {
            rosenbluth_sample_helix(region_size, o_h, u, v, excluded_volume, helix_R_factor);
            generate(region_size, v, u, o_h, s_x, s_y, rcentres);

        }
        else {
            sample_helix_vectors(u, v);
            generate(region_size, v, u, o_h, s_x, s_y, rcentres);
            // when rosenbluth sampling is off, the excluded volume interaction between monomers that are not regrown and the monomers 
            // that form the helix must be accounted for in the acceptance rule. (otherwise the term cancels).
            std::vector<std::vector<double>> helix_pos{ s_x };
            helix_pos.insert(helix_pos.end(), s_y.begin(), s_y.end());
            helix_excluded_volume_interaction(helix_pos);

        }

        // grow the new configuration (the linked one)
        int gamma{ static_cast<int>(pow(-1.0,alpha + beta)) };
        for (int i{ 0 }; i < region_size; i++) {
            old_config_positions[origin + gamma * i] = s_x[i];
            old_config_positions[origin_y + gamma * i] = s_y[i];

            excluded_volume.push_back(s_x[i]);
            excluded_volume.push_back(s_y[i]);
        }


        grow_limits(link_limits, alpha, false);
        old_growth_lims;
        //need to do helix excluded volume 
        print_1d_int_vec(destructure);

    }
    //trying to catch a bug where there's a huge separation
    for (auto i : chain) {
        int j = i->get_id();
        if (j + 1 > chain.size() - 1) {
            break;
        }
        std::vector<double> temp1{ chain[j + 1]->get_position() }, temp2{ i->get_position() };
        if (dist_2_points3d(temp1, temp2) > 5) {
            std::cout << "caught one" << std::endl;
        }
    }

}

void polymer::unlink_update(int s_i) {

    std::vector<int> s{ helix_list[s_i]->get_monomers() };
    //add structure that we unlinked back to the pool of possible regions for a link move
    search_results.push_back(s);

    // update the unlinked monomers. they are not part of a structure anymore. structure member variable == 0 if not linked.
    unpack_region(s);
    for (auto i : s) {
        chain[i]->change_structure(0);
    }

    //delete corresponding helix object and delete element from array.
    delete helix_list[s_i];
    helix_list.erase(helix_list.begin() + s_i);

    // update positions and neighbours
    update_positions();

}

void polymer::unlink_growth_limits(std::vector<int> helix, int alpha, int beta, std::vector<std::vector<int>>& growth_limits)
{
    std::vector<int> reaction_region{ helix };
    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end(), std::greater<>());

    }

    int a, b, c, d;

    a = reaction_region[0];
    b = reaction_region[1];
    c = reaction_region[2];
    d = reaction_region[3];

    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end());

    }

    //define w,x,y,z as the leftmost linker, middle linkers and rightmost linkers respectively in the most
    //complicated case where there are structures on either side of the structure to be unlinked as well as in
    //the middle.
    int w, x, y, z;//sometimes some of these will be the same.
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];


    int region_size{ abs(b - a) + 1 };
    bool middle_struct{ true };
    if (abs(b - c) < abs(b - x)) {
        //x = c;
        y = b;
        middle_struct = false;
    }

    std::vector<std::vector<int>> degrowth_limits;

    if (beta == 0) {
        degrowth_limits.push_back({ a,x });
        if (middle_struct) {
            degrowth_limits.push_back({ y,z });
        }
    }
    else if (beta == 1) {
        if (middle_struct) {
            degrowth_limits.push_back({ w,b });
            degrowth_limits.push_back({ y,z });
        }
        else {
            degrowth_limits.push_back({ w,z });
        }

    }
    growth_limits = degrowth_limits;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool polymer::linked()
{
    if (helix_list.size() != 0) {
        return true;
    }
    else return false;
}

void polymer::link(std::vector<int>& link_region, int alpha, int beta, int ss_index)
{
    std::cout << "alpha beta " << alpha << " " << beta << std::endl;
    //neighbouring_linkers();//update neighbours
    std::vector<int> reaction_region{ link_region };



    int kappa{ static_cast<int>(pow(-1.0,alpha)) };
    int origin{ reaction_region[alpha * 3 + kappa * beta] }, origin_y{ reaction_region[(1 - alpha) * 3 - kappa * (1 - beta)] };

    // define a,b,c,d as the extremes of the structure to be destroyed. either a,b,c,d could be the origin of the unlink.
    int a, b, c, d;
    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end(), std::greater<>());

    }

    a = reaction_region[0];
    b = reaction_region[1];
    c = reaction_region[2];
    d = reaction_region[3];

    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end());

    }

    //define w,x,y,z as the leftmost linker, middle linkers and rightmost linkers respectively in the most
    //complicated case where there are structures on either side of the structure to be unlinked as well as in
    //the middle.
    int w, x, y, z;//sometimes some of these will be the same.
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    int region_size{ abs(b - a) + 1 };

    if (abs(b - c) <= abs(b - x)) {
        x = c;
        y = b;
    }
    std::vector<std::vector<int>> link_limits;

    if (beta == 1) {
        link_limits.push_back({ w,a });
    }
    else if (beta == 0 && y != b) {
        link_limits.push_back({ b,x });
    }
    link_limits.push_back({ y,c });
    link_limits.push_back({ d,z });

    print_1d_int_vec(reaction_region);
    std::cout << " link growth limits" << std::endl;
    print_2d_int_vec(link_limits);


    //start the growth. 

    //first update the excluded volume
    update_excluded_volume(link_limits, reaction_region);
    std::cout << "excluded volume size at very start " << excluded_volume.size() << std::endl;
    // grow helix

    std::vector<std::vector<double>> s_x, s_y;
    std::vector<double> u, v;
    std::vector<double> o_h{ chain[origin]->get_position() };
    std::vector<std::vector<double>> rcentres(2);


    // if we've turned on rosenbluth sampling then we rosenbluth sample the helix. the rosenbluth sample helix
    // returns the selected helix vectors by reference.
    if (rosenbluth_switch == true) {
        rosenbluth_sample_helix(region_size, o_h, u, v, excluded_volume, helix_R_factor);
        generate(region_size, v, u, o_h, s_x, s_y, rcentres);

    }
    else {
        sample_helix_vectors(u, v);
        generate(region_size, v, u, o_h, s_x, s_y, rcentres);
        // when rosenbluth sampling is off, the excluded volume interaction between monomers that are not regrown and the monomers 
        // that form the helix must be accounted for in the acceptance rule. (otherwise the term cancels).
        std::vector<std::vector<double>> helix_pos{ s_x };
        helix_pos.insert(helix_pos.end(), s_y.begin(), s_y.end());
        helix_excluded_volume_interaction(helix_pos);

    }

    // grow the new configuration (the linked one)
    int gamma{ static_cast<int>(pow(-1.0,alpha + beta)) };
    for (int i{ 0 }; i < region_size; i++) {
        new_config_positions[origin + gamma * i] = s_x[i];
        new_config_positions[origin_y + gamma * i] = s_y[i];

        excluded_volume.push_back(s_x[i]);
        excluded_volume.push_back(s_y[i]);
    }
    grow_limits(link_limits, alpha,true);
    proposed_link_helix = new helix_struct(reaction_region, v, rcentres, alpha, beta);

    // grow the old configuration 
    std::vector<std::vector<int>> unlink_limits;
    unlink_growth_limits(reaction_region, alpha, beta, unlink_limits);
    update_excluded_volume(unlink_limits);
    std::cout << "unlink growth limits" << std::endl;
    print_2d_int_vec(unlink_limits);
    grow_limits(unlink_limits, alpha, false);

}

void polymer::link_update(int ss_index) {
    //, int a, int b, std::vector<double>&v_vec, std::vector<std::vector<double>> &rcentres
    std::vector<int> s{ search_results[ss_index] };

    helix_list.push_back(proposed_link_helix);

    unpack_region(s);
    for (auto i : s) {
        chain[i]->change_structure(helix_list.size() + 1);
    }
    search_results.erase(search_results.begin() + ss_index);//this line ensures that we remove the linked structure from the pool of

    //neighbouring_linkers();//update neighbours


    update_positions();

}

void polymer::link_growth_limits(std::vector<int> helix, int alpha, int beta, std::vector<std::vector<int>>& growth_limits)
{
    std::vector<int> reaction_region{ helix };

    // define a,b,c,d as the extremes of the structure to be destroyed. either a,b,c,d could be the origin of the unlink.
    int a, b, c, d;
    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end(), std::greater<>());

    }

    a = reaction_region[0];
    b = reaction_region[1];
    c = reaction_region[2];
    d = reaction_region[3];

    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end());

    }

    //define w,x,y,z as the leftmost linker, middle linkers and rightmost linkers respectively in the most
    //complicated case where there are structures on either side of the structure to be unlinked as well as in
    //the middle.
    int w, x, y, z;//sometimes some of these will be the same.
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    int region_size{ abs(b - a) + 1 };

    if (abs(b - c) <= abs(b - x)) {
        x = c;
        y = b;
    }
    std::vector<std::vector<int>> free_growth_limits;

    if (beta == 1) {
        free_growth_limits.push_back({ w,a });
    }
    else if (beta == 0 && y != b) {
        free_growth_limits.push_back({ b,x });
    }
    free_growth_limits.push_back({ y,c });
    free_growth_limits.push_back({ d,z });

    growth_limits = free_growth_limits;
}


bool polymer::reject_link(std::vector<int>& link_region, int alpha, int beta)
{
    if (link_region.empty()) {
        return true;
    }
    if (overlapping_monomers(link_region)) {
        return true;
    }
    if (helix_list.size() == 0) {
        return false;
    }
    std::vector<int> reaction_region{ link_region };
    // s(ab) = [00,01,11,10] [0, 1, 2 ,3]
    int kappa{ static_cast<int>(pow(-1.0,alpha)) };
    int origin{ reaction_region[alpha * 3 + kappa * beta] };

    // define a,b,c,d as the extremes of the structure to be destroyed. either a,b,c,d could be the origin of the unlink.
    int a, b, c, d;
    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end(), std::greater<>());

    }

    a = reaction_region[0];
    b = reaction_region[1];
    c = reaction_region[2];
    d = reaction_region[3];

    if (alpha == 1) {
        std::sort(reaction_region.begin(), reaction_region.end());

    }

    //define w,x,y,z as the leftmost linker, middle linkers and rightmost linkers respectively in the most
    //complicated case where there are structures on either side of the structure to be unlinked as well as in
    //the middle.
    int w, x, y, z;//sometimes some of these will be the same.
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    //std::cout << "abcd " << a << " " << b << " " << c << " " << d << std::endl;
    //std::cout << "wxyz " << w << " " << x << " " << y << " " << z << std::endl;

    int region_size{ abs(b - a) + 1 };
    bool middle_struct{ true };
    if (abs(b - c) < abs(b - x)) {// this means that there's no structure in between
        x = c;
        y = b;
        middle_struct = false;
    }

    // we always accept the first structure.
    if (helix_list.size() == 0) {
        return false;
    }

    std::vector<double> temp1{ chain[origin]->get_position() }, temp2(3);
    if (beta == 1) {
        if (z != (1 - alpha) * (chain.size() - 1)) {// there is a structure on the extreme end
            temp1 = chain[origin]->get_position();
            temp2 = chain[z]->get_position();
            if (dist_2_points3d(temp1, temp2) + sideways_length(region_size)
    > abs(d - z)) {
                return true;

            }
        }

        if (w != alpha * (chain.size() - 1)) {
            temp2 = chain[w]->get_position();
            if (dist_2_points3d(temp1, temp2)
                + side_length(region_size) > abs(w - a)) {
                //std::cout << "constraint 1 violated" << std::endl;

                return true;
            }
        }

        if (abs(b - c) > abs(b - x)) {// there is a structure in between
            temp2 = chain[y]->get_position();
            if (dist_2_points3d(temp1, temp2) + helix_separation()
        > abs(y - c)) {
                std::cout << "constraint check with helix separation: " << y << " " << c << std::endl;
                return true;
            }

        }
        return false;
    }
    else if (beta == 0) {

        if (z != (1 - alpha) * (chain.size() - 1)) {

            temp2 = chain[z]->get_position();
            if (dist_2_points3d(temp1, temp2) + helix_separation()
                > abs(d - z)) {
                return true;
            }

        }

        if (abs(b - c) > abs(b - x)) {

            temp2 = chain[x]->get_position();

            if (dist_2_points3d(temp1, temp2)
                + side_length(region_size) > abs(b - x)) {
                //std::cout << "constraint 1 violated" << std::endl;

                return true;
            }

            temp2 = chain[y]->get_position();
            if (dist_2_points3d(temp1, temp2) + sideways_length(region_size)
            > abs(y - c)) {
                return true;

            }

        }
        return false;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double polymer::link_acceptance(bool link_move)
{
    // boolean argument link_move
    // link_move == 0 if we are doing an unlink move
    // link_move == 1 if we are doing a link move

    double W_b, W_u, Z_ni_A{ 0.27 }, Z_ni_B{ 0.27 }, Z_ni_AB{ 1 }, G0{ -0.5}, rho0{ 1 }, U_b, U_u;
    double beta{ 1 };


     // selection probability. for a link move we have to select a set of compatible regions randomly
    // for an unlink we have to select an existing double helix structure randomly.
    // call this P_preselection to distinguish it from the selection due to Rosenbluth sampling.
    double P_preselection, Po2n, Pn2o;// n2o = new to old. o2n vice versa.

    if (link_move) {// link is the proposed move
        Po2n = 1 / static_cast<double>(search_results.size()), Pn2o = 1 / (static_cast<double>(helix_list.size()) + 1);

    }
    else if (!link_move) {// unlink is the proposed move
        Po2n = 1 / static_cast<double>(helix_list.size()), Pn2o = 1 / (static_cast<double>(search_results.size()) + 1);
    }
    P_preselection = Pn2o / Po2n;


    double z{ 1 };
    double W_gen_old, W_gen_new, W_helix{helix_R_factor};

    // the acceptance rules only feature the weights. 
    W_gen_old = weight_gen_config(false);
    W_gen_new = weight_gen_config(true);

    // if P_gen_new is zero because a monomer was grown with weight 0 (actually we would get a division by 0 but we explicitly tell p_gen_configuration function to return 0 in the case 
    // of a zero weight monomer) we reject the move (acceptance = 0).
    if (W_gen_new == 0.0) {
        return 0;
    }
    else {
        double prefactors{ (exp(-beta * G0) / rho0) * (Z_ni_A * Z_ni_B / Z_ni_AB) };
        if (!link_move) {
            prefactors = 1 / prefactors;
            W_helix = 1 / W_helix;

        }
        old_weights;
        new_weights;
        old_growth_lims;
        new_growth_lims;
        double W_gen_ratio{ W_gen_new* W_helix / W_gen_old };
        z = W_gen_ratio * prefactors * P_preselection;
        if (std::isnan(z) || std::isinf(z)) {
            std::cout << "stop" << std::endl;

        }
        std::cout << " link acceptance ratio z " << z << std::endl;
        if (z >= 1) {
            return 1;
        }
        else {
            return z;
        }


    }

}



double polymer::zip_acceptance(bool zip_move, int side)
{
    // if zip_move == 1 then we are doing a zip move
    // if zip_move == 0 then we are doing an unzip move


    // the probability of picking a zipped structures to be unzipped
    double P_uz{ 1 / static_cast<double>(zipped_structures.size()) };
    // the probability of picking a structure to be zipped from the set of zippable structures
    double P_zs{ 1 / static_cast<double>(extendable_structures.size()) };

    double Z1_Mplus1{ 0.5 }, Z1_M{ 0.5 }, Z2_Mplus1{ 0.5 }, Z2_M{ 0.5 };
    double dG_M{ 1 }, dG_Mplus1{ 1.1 };
    double deltadG{ -(dG_Mplus1 - dG_M) };

    if (new_growth_lims.size() > 1 && (new_growth_lims[0][0] > new_growth_lims[1][0])) {
        std::reverse(new_growth_lims.begin(), new_growth_lims.end());
    }
    
    std::vector<std::vector<int>> old_growth_lims(new_growth_lims.size());
    int kappa;
    zip_move == 1 ? kappa = -1 : kappa = +1;
    if (new_growth_lims.size() == 1) {
        int i1{ new_growth_lims[0][0] }, i2{ new_growth_lims[0][1] };
        old_growth_lims[0] = { i1 + kappa, i2 - kappa };
    }
    else {
        std::vector<int> temp(2);
        if (side == 0) {
            temp[0] = new_growth_lims[0][0];
            temp[1] = new_growth_lims[0][1] - kappa;
            old_growth_lims[0] = temp;
            //old_growth_lims[0][0] = regrowth_lims[0][0];
            //old_growth_lims[0][1] = regrowth_lims[0][1] - kappa;

            temp[0] = new_growth_lims[1][0] + kappa;
            temp[1] = new_growth_lims[1][1];
            old_growth_lims[1] = temp;

            //old_growth_lims[1][0] = regrowth_lims[1][0] + kappa;
            //old_growth_lims[1][1] = regrowth_lims[1][1];
        }
        else if (side == 1) {
            temp[0] = new_growth_lims[0][0]+kappa;
            temp[1] = new_growth_lims[0][1];
            old_growth_lims[0] = temp;

            //old_growth_lims[0][0] = regrowth_lims[0][0] + kappa;
            //old_growth_lims[0][1] = regrowth_lims[0][1];

            temp[0] = new_growth_lims[1][0];
            temp[1] = new_growth_lims[1][1] - kappa;
            old_growth_lims[1] = temp;

            //old_growth_lims[1][0] = regrowth_lims[1][0];
            //old_growth_lims[1][1] = regrowth_lims[1][1] - kappa;

        }
    }    
    std::vector<double> unzipped_probs(new_growth_lims.size()), zipped_probs(new_growth_lims.size());

    // in this block we get the rosenbluth weights that we need for the sections of interest
    for (int i = 0; i < new_growth_lims.size(); i++)
    {
        if (zip_move == 0) {
            zipped_probs[i] = old_weights->subsection_probability(old_growth_lims[i]);
            unzipped_probs[i] = new_weights->subsection_probability(new_growth_lims[i]);

        }
        else {
            zipped_probs[i] = new_weights->subsection_probability(new_growth_lims[i]);
            unzipped_probs[i] = old_weights->subsection_probability(old_growth_lims[i]);

        }

    }

    // in this block we get the fixed endpoint probability factor if the limit of interest is between
    // 2 fixed endpoints.
    std::vector<int> limit(2);
    for (int i{ 0 }; i < new_growth_lims.size(); i++) {
        limit = new_growth_lims[i];
        if (limit[0] != 0 && limit[1] != chain.size() - 1) {
            std::vector<double> temp1(3), temp2(3);
            temp1 = chain[limit[0]]->get_position();
            temp2 = chain[limit[1]]->get_position();

            temp1 = vector_subtraction(temp1, temp2);
            double temp{ vector_modulus(temp1) };

            if (abs(limit[1] - limit[0]) < 2) {
                //unzipped_probs[i] *= 1 / 2 * pi;
                //zipped_probs[i] *= 1 / 2 * pi;
                unzipped_probs[i] = unzipped_probs[i] / (8 * atan(1));

                zipped_probs[i] = zipped_probs[i] / (8 * atan(1));

            }
            else {

                if (!zip_move) {
                    unzipped_probs[i] = unzipped_probs[i] * ideal_chain_pdf(temp, abs(limit[1] - limit[0]));

                }
                else if (zip_move) {
                    zipped_probs[i] = zipped_probs[i] * ideal_chain_pdf(temp, abs(limit[1] - limit[0]));
                }

            }
        }
    }
    for (int i{ 0 }; i < old_growth_lims.size(); i++) {
        limit = old_growth_lims[i];
        if (limit[0] != 0 && limit[1] != chain.size() - 1) {
            std::vector<double> temp1(3), temp2(3);
            temp1 = chain[limit[0]]->get_position();
            temp2 = chain[limit[1]]->get_position();

            temp1 = vector_subtraction(temp1, temp2);
            //check if this is just helix separation.
            double temp{ vector_modulus(temp1) };


            if (abs(limit[1] - limit[0]) < 2) {
                if (!zip_move) {
                    zipped_probs[i] = zipped_probs[i] / (8 * atan(1));

                }
                else {
                    unzipped_probs[i] = unzipped_probs[i] / (8 * atan(1));
                }


            }
            else {
                if (!zip_move) {
                    zipped_probs[i] = zipped_probs[i] * ideal_chain_pdf(temp, abs(limit[1] - limit[0]));

                }
                else if (zip_move) {
                    unzipped_probs[i] = unzipped_probs[i] * ideal_chain_pdf(temp, abs(limit[1] - limit[0]));
                }
            }

        }
    }

    double P_gen_unzipped{product_of_elements(unzipped_probs)}, P_gen_zipped{product_of_elements(zipped_probs)};

    double z;

    //z = exp(-deltadG) * (pow(Z1_Mplus1, 2) / pow(Z1_M, 2)) * (Z2_M / Z2_Mplus1) * (P_gen_unzipped / P_gen_zipped);
    z = exp(-1.0*deltadG) * (pow(Z1_Mplus1, 2) / pow(Z1_M, 2)) * (Z2_M / Z2_Mplus1) * (P_gen_zipped / P_gen_unzipped);
    if (std::isinf(z)) {
        std::cout << "z infinity" << std::endl;
    }
    if (!zip_move) { z = 1 / z; }
    std::cout << "zip acceptance " << z << std::endl;
    if (z != z) {
        std::cout << "z is undefined" << std::endl;
    }

    if (z >= 1) {
        return 1;
    }
    else {
        return z;
    }
}

double polymer::swivel_acceptance()
{
    double z{1};
    for (auto i : new_growth_lims) {
        z = z * p_gen_ratio(i);
        if (std::isinf(z)) {
            std::cout << "caught" << std::endl;
        }
    }
    //double z = old_config_prob / new_config_prob;
    
    if (std::isinf(z)) {
        std::cout << "z infinity" << std::endl;
    }

    std::cout << "swivel acceptance ratio " << z << std::endl;
    if (z >= 1) {
        return 1;
    }
    else { return z; }

}


double polymer::weight_gen_config(bool forward_move)
{

    double weight_gen{ 1 };
    std::vector<int> limit(2), regrown_indices(2);
    std::vector<double> rA(3), rB(3), rAB(3);
    double mod_rAB;

    std::vector<std::vector<int>>* limits;
    std::vector<std::vector<double>>* positions;
    rosenbluth_growth* weights;

    if (forward_move == true) {
        weights = new_weights;
        limits = &new_growth_lims;
        positions = &new_config_positions;

        // if one of the monomers of the new configuration was generated with weight 0 then we need to reject the move. if we dont return 0 here we'll get a nan. 
        // another solution would be to always immediately calculate the ratio of old to new.
        for (auto i : weights->get_weights()) {
            if (i == 0.0) {
                return 0.0;
            }
        }
    }
    else {
        weights = old_weights;
        limits = &old_growth_lims;
        positions = &old_config_positions;

    }

    for (auto i : *limits) {
        limit = i;
        if (limit[0] > limit[1]) {
            std::reverse(limit.begin(), limit.end());
        }
        int A{ limit[0] }, B{ limit[1] };

        if (abs(B - A) == 1) {
            continue;
        }

        if (limit[0] != 0 && limit[1] != chain.size() - 1) {// for both ends fixed. relevant for swivel and link.
            // will give the ratio of new to old
            regrown_indices = { A + 1,B - 1 };// we don't take the weights of the endpoints

            weight_gen *= weights->subsection_weight(regrown_indices);// rosenbluth distribution probability for section considered.

            int n_segments{ abs(B - A) };

            rA = positions->operator[](A), rB = positions->operator[](B);
            rAB = vector_subtraction(rA, rB);
            mod_rAB = vector_modulus(rAB);

            weight_gen *= ideal_chain_pdf(mod_rAB, n_segments); // "yamakawa" probability

        }
        else { // random walk branch
            if (A == 0) {
                regrown_indices = { A,B - 1 };
            }
            else if (B == chain.size() - 1) {
                regrown_indices = { A + 1,B };
            }
            weight_gen *= weights->subsection_weight(regrown_indices);// rosenbluth distribution probability for section considered.
        }

    }

    return weight_gen;
}

double polymer::p_gen_ratio(std::vector<int> limit)// the ratio of p(new to old) / p(old to new) for the generation of new subsections.
// assumes that the growth limits are the same in both cases. this is true for link, unlink and swivel  (but not for zip/unzip??)
{
    double ratio;
    double R_old, R_new;
    std::vector<int> regrown_indices(2);
    if (limit[0] > limit[1]) {
        std::reverse(limit.begin(), limit.end());
    }
    int A{ limit[0] }, B{ limit[1] };

    if (abs(B - A) == 1) {
        return 1.0;
    }

    if (limit[0] != 0 && limit[1] != chain.size() - 1) {// for both ends fixed. relevant for swivel and link.
        // will give the ratio of new to old
        // 


        int n_segments{ abs(B - A) };
        regrown_indices = { A + 1,B - 1 };// we don't take the weights of the endpoints
        R_old = old_weights->subsection_probability(regrown_indices);// should calculate this differently in case R_new has a 0 weight in it.
        // in that case there will be division by 0. 
        R_new = new_weights->subsection_probability(regrown_indices);

        if (R_new == 0.0) {
            return 0.0;
        }

        std::vector<double> roA{ chain[A]->get_position() }, roB{ chain[B]->get_position() },
            rnA{ new_config_positions[A] }, rnB{ new_config_positions[B] },
            roE2E{ vector_subtraction(roA,roB) }, rnE2E{ vector_subtraction(rnA,rnB) };

        double roAB{ vector_modulus(roE2E) }, rnAB{ vector_modulus(rnE2E) };

        ratio = (R_new / R_old) * (ideal_chain_pdf(rnAB, n_segments)/ideal_chain_pdf(roAB, n_segments));

    }
    else {
        if (A == 0) {
            regrown_indices = { A,B - 1 };
        }
        else if (B == chain.size() - 1) {
            regrown_indices = { A + 1,B };
        }
        R_old = old_weights->subsection_probability(regrown_indices);
        R_new = new_weights->subsection_probability(regrown_indices);
        if (R_new == 0) {
            return 0.0;
        }

        ratio = (R_new / R_old);
    }

    return ratio;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void polymer::structure_extension(std::vector<int>& s, int& index) {

    // label s_3 and s_5 the 3' and 5' compatible regions forming the structure s
    // s_3 and s_5 themselves have 3' and 5' ends

    // search compatibility between s_3 3' and s_5 5' ends
    int s1{ s.front() }, s2{ s.back() };
    int lim;
    (chain.size() - 1) - s2 > s1 ? lim = s1 : lim = (chain.size() - 1) - s2;

    std::string b1, b2;

    std::vector<int> c3, c5;
    std::vector<int> extended_region, temp1;
    int i1, i2;
    for (int i{ 1 }; i < lim; i++) {
        i1 = s1 - i, i2 = s2 + i;
        b1 = chain[s1 - i]->get_base();
        b2 = chain[s2 + i]->get_base();
        if (regions_compatible(b1, b2) != true) {
            break;
        }
        else {
            //std::cout << "it can be extended1 !" << std::endl;
            c3.push_back(s1 - i);
            c5.push_back(s2 + i);
        }
    }
    if (c3.size() > 0) {
        temp1.push_back(c3.back());
        temp1.push_back(c5.back());

    }
    else {
        temp1.push_back(s.front());
        temp1.push_back(s.back());
    }

    c3.insert(c3.end(), c5.begin(), c5.end());
    extended_region = c3;
    // search compatibility between s_3 5' and s_5 3' ends
    c3.clear(), c5.clear();
    s1 = s[s.size() / 2 - 1], s2 = s[s.size() / 2];
    lim = floor((s2 - s1) / 2);


    for (int i{ 1 }; i < lim; i++) {
        i1 = s1 + i, i2 = s2 - i;

        b1 = chain[s1 + i]->get_base();
        b2 = chain[s2 - i]->get_base();
        if (regions_compatible(b1, b2) != true) {
            break;
        }
        else {
            //std::cout << "it can be extended2 !" << std::endl;
            c3.push_back(s2 - i);
            c5.push_back(s1 + i);
        }
    }
    if (c3.size() > 0) {
        temp1.push_back(c3.back());
        temp1.push_back(c5.back());

    }
    else {
        temp1.push_back(s1);
        temp1.push_back(s2);
    }

    c3.insert(c3.end(), c5.begin(), c5.end());
    extended_region.insert(extended_region.end(), c3.begin(), c3.end());
    extended_region.insert(extended_region.end(), s.begin(), s.end());
    std::sort(extended_region.begin(), extended_region.end());

    std::sort(temp1.begin(), temp1.end());
    pack_region(extended_region);

    int c{ 0 };
    for (auto i : extended_region) {
        for (auto j : s) {
            if (i == j) {
                c++;
            }
        }
    }
    if (c == 0) {
        std::cout << "extendable on both sides" << std::endl;
    }
    s = extended_region;
    //print_1d_int_vec(extended_region);
    //print_1d_int_vec(temp1);
    //print_bases1();
}

bool polymer::zip_structure_overlap(std::vector<int> s_z, int side) {
    if (side == 0) {
        if (chain[s_z[0]]->get_structure() != 0) {
            return true;
        }
        if (chain[s_z[3]]->get_structure() != 0) {
            return true;
        }


    }
    else if (side == 1) {
        if (chain[s_z[1]]->get_structure() != 0) {
            return true;
        }
        if (chain[s_z[2]]->get_structure() != 0) {
            return true;
        }

    }
    return false;
}


bool polymer::reject_zip(std::vector<int>& extension, std::vector<std::vector<double>> new_positions, int side) {
    int w, z;

    if (side == 0) {
        w = chain[extension[0]]->left_right_link[0];
        z = chain[extension[3]]->left_right_link[1];

        if (w != 0) {
            std::vector<double> temp{ chain[w]->get_position() };
            if (dist_2_points3d(new_positions[0], temp)>abs(w-extension[0])) {
                return true;
            }
        }

        if (z != chain.size() - 1) {
            std::vector<double> temp{ chain[z]->get_position() };
            if (dist_2_points3d(new_positions[1], temp) > abs(z - extension[3])) {
                return true;
            }
        }
    }
    else if (side == 1) {
        w = chain[extension[1]]->left_right_link[1];
        z = chain[extension[2]]->left_right_link[0];
        if (w < extension[2]) {
            std::vector<double> temp{ chain[w]->get_position() };
            if (dist_2_points3d(new_positions[0], temp) > abs(extension[1] - w)) {
                return true;
            }
        }
        if (z > extension[1]) {
            std::vector<double> temp{ chain[z]->get_position() };
            if (dist_2_points3d(new_positions[1], temp) > abs(z - extension[2])) {
                return true;
            }
        }            


    }
    return false;
}


void polymer::zip_growth_limits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits) {

    std::vector<int> reaction_region{ s };
    std::vector<std::vector<int>> limits;

    // define a,b,c,d as the extremes of the structure to be destroyed. either a,b,c,d could be the origin of the unlink.
    int a, b, c, d;
    a = reaction_region[0];
    b = reaction_region[1];
    c = reaction_region[2];
    d = reaction_region[3];

    int p, q;
    if (side == 0) {
        p = a;
        q = d;
    }
    else {
        p = b;
        q = c;
    }


    //define w,x,y,z as the leftmost linker, middle linkers and rightmost linkers respectively in the most
    //complicated case where there are structures on either side of the structure to be unlinked as well as in
    //the middle.
    int w, x, y, z;//sometimes some of these will be the same.
    w = chain[a]->left_right_link[0];
    x = chain[b]->left_right_link[1];
    y = chain[c]->left_right_link[0];
    z = chain[d]->left_right_link[1];

    int region_size{ abs(b - a) + 1 };
    bool middle_struct{ true };
    if (abs(b - c) <= abs(b - x)) {// this means that there's no structure in between
        x = c;
        y = b;
        middle_struct = false;
    }

    if (side == 0) {
        limits.push_back({ w,p });
        limits.push_back({ q,z });
    }
    else if (side == 1) {
        limits.push_back({ p,x });
        if (middle_struct) {
            limits.push_back({ y,q });
        }
    }

    growth_limits = limits;
}

void polymer::unzip_growth_linits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits) {
    std::vector<int> reaction_region{ s };
    std::vector<std::vector<int>> limits;
    int a, b, c, d;
    a = reaction_region[0];
    b = reaction_region[1];
    c = reaction_region[2];
    d = reaction_region[3];

    if (side == 0) {
        int w, z;
        w = chain[a]->left_right_link[0];
        z = chain[d]->left_right_link[1];

        int a1{ a + 1 }, d1{ d - 1 };
        limits.push_back({ w,a1 });
        limits.push_back({ d1,z });
    }
    else if (side == 1) {
        int x{ chain[b]->left_right_link[1] };
        int b1{ b - 1 }, c1{ c + 1 };
        if (abs(b - x) < abs(b - c)) {
            int y{ chain[c]->left_right_link[0] };
            limits.push_back({ b1,x });
            limits.push_back({ y,c1 });

        }
        else {
            limits.push_back({ b1,c1 });
        }
    }
    growth_limits = limits;
}


void polymer::sample_zip_move(bool &success, int&s_index) {
    //neighbouring_linkers();//update neighbours

    //cannot perform zip if none of the linked structures can be extended.
    update_extensible_structures();

    if (extendable_structures.size() == 0) {
        success = false;
        return;
    }

    int structure_index{ static_cast<int>(rand2(0,extendable_structures.size())) };
    s_index = structure_index;

}
void polymer::zip(bool& success, int &sigma, int &s_index) {

    //neighbouring_linkers();//update neighbours

    //cannot perform zip if none of the linked structures can be extended.
    update_extensible_structures();
    //std::cout << "got out of update_extensible structures" << std::endl;
    if (extendable_structures.size() == 0) { 
        success = false;
        return; }


    // sample a structure from extendable structures. p_select = 1 / extendable_structures.size()
    int structure_index{ static_cast<int>(rand2(0,extendable_structures.size())) };
    structure_index = extendable_structures[structure_index];
    helix_struct* extend_helix{ helix_list[structure_index] };
    std::vector<int> small_helix{ extend_helix->get_monomers() }, s_z{ small_helix };
    std::vector<int> reaction_region = extend_helix->get_extension();
    //sample side
    int side{ static_cast<int>(rand2(0,2)) };
    sigma = side;
    s_index = structure_index;

    sample_zip_region(s_z, reaction_region,side);
    print_1d_int_vec(s_z);

    if (s_z.size() == 0) {//cant zip a double helix if there are none
        success = 0;
        return;
    }
    if (s_z[0]<0 || s_z[3]>chain.size() - 1) { // can't try to zip a region which is outside of the terminal endpoints of the polymer
        std::cout << "invalid zip region" << std::endl;
        success = 0;
        return;
    }
    if (zip_structure_overlap(s_z, side)) { // we don't allow zipping structures which would overlap with other structures (ie a merging of structures).
        success = 0;
        return;
    }
    std::vector<std::vector<double>> new_pos(2);
    new_pos = add_bp_to_helix(extend_helix, side);
    print_2d_doub_vec(new_pos);
    if (reject_zip(s_z, new_pos, side)) {
        std::cout << "zip move rejected due to connectivity constraint or structure overlap." << std::endl;
        success = 0;
        return;
    }
    print_1d_int_vec(s_z);

    std::vector<std::vector<int>> growth_limits;
    zip_growth_limits(s_z, side, growth_limits);
    update_excluded_volume(growth_limits, s_z);

    //have to repeat some of the code from add bp to helix now that the zip has been accepted.
    // zip has been accepted here based on physical constraints but not yet by the acceptance probability.
    // this section will need to be moved to a separate function so that we only update positions and other
    // data when a move has been accepted.
    //********************************************************
    int o_x, o_y;

    if (side == 0) {
        o_x = small_helix[0];
        o_y = small_helix[3];
    }
    else if (side == 1) {
        o_x = small_helix[1];
        o_y = small_helix[2];

    }
    std::vector<double> running_centre(3);
    int beta{ extend_helix->get_beta() };

    side == 0 ? running_centre = extend_helix->get_rc(beta) : running_centre = extend_helix->get_rc(1 - beta);
    double pitch{ 0.28 };
    if (abs((1 - beta) - (1 - side)) == 0) {
        pitch = -pitch;
    }
    std::vector<double> v{ extend_helix->get_v() };

    std::vector<double> temp(3);
    temp = multiplication_by_scalar(pitch, v);
    running_centre = vector_addition(running_centre, temp);

    //side == 0 ? extend_helix->set_rc(beta, running_centre) : extend_helix->set_rc(1 - beta, running_centre);
    
    new_running_centre = running_centre;
    zip_unzip_structure = s_z;
    //********************************************************

    // change the positions of the zipped monomers
    int kappa{ static_cast<int>(pow(-1,side)) };
    new_config_positions[o_x - kappa] = new_pos[0];
    new_config_positions[o_y + kappa] = new_pos[1];


    // have to add all the monomers in the helix back to the excluded volume. also have to calculate the interaction of the helix with
    // the ungrown (previous) sections.
    std::vector < std::vector<double>> helix_positions(2);
    for (auto i : small_helix) {
        excluded_volume.push_back(chain[i]->get_position());
    }
    helix_positions[0] = new_pos[0];
    helix_positions[1] = new_pos[1];

    helix_excluded_volume_interaction(helix_positions);

    excluded_volume.push_back(new_pos[0]);
    excluded_volume.push_back(new_pos[1]);

    //extend_helix->extend(s_z);

    grow_limits(growth_limits, 0,true);// this function updates the Rosenbluth weights 

    //trying to catch a bug where there's a huge separation
    for (auto i : chain) {
        int j = i->get_id();
        if (j + 1 > chain.size() - 1) {
            break;
        }
        std::vector<double> temp1{ chain[j + 1]->get_position() }, temp2{ i->get_position() };
        if (dist_2_points3d(temp1, temp2) > 5) {
            std::cout << "caught one" << std::endl;
            std::cout << "id " << i->get_id() << std::endl;
        }
    }
    std::cout << "sigma " << side << std::endl;

    std::cout << "ZIP MOVE SUCCESS" << std::endl;
    success = true;

}
void polymer::zip_update(int s_index, int side) {
    int beta{ helix_list[s_index]->get_beta() };
    side == 0 ? helix_list[s_index]->set_rc(beta, new_running_centre) : helix_list[s_index]->set_rc(1 - beta, new_running_centre);
    helix_list[s_index]->extend(zip_unzip_structure);
    update_positions();

    unpack_region(zip_unzip_structure);
    for (auto i : zip_unzip_structure) {
        chain[i]->change_structure(s_index + 1);
    }


}
void polymer::unzip(bool &success, int &sigma, int &s_index) {
    //neighbouring_linkers();//update neighbours

    update_large_struct_list();
    if (zipped_structures.size() == 0) {
        //std::cout << "unzip move not possible, no regions to unzip" << std::endl;
        success = 0;
        return;
    }


    // sample a structure from the extended structures. p_select = 1 / extended_structures.size()
    int structure_index{ static_cast<int>(rand2(0,zipped_structures.size())) };
    structure_index = zipped_structures[structure_index];
    std::vector<int> unzip_this{ helix_list[structure_index]->get_monomers() };
    std::vector<int> zipped_helix{ unzip_this };

    
    //sample side
    int side{ static_cast<int>(rand2(0,2)) };
    sample_unzip_region(unzip_this, side);
    //std::cout << "unzip went through: region" << std::endl;
    print_1d_int_vec(unzip_this);
    std::vector<std::vector<int>> limits;
    unzip_growth_linits(zipped_helix, side, limits);
    sigma = side;

    update_excluded_volume(limits);
    grow_limits(limits, 0,true);

    zip_unzip_structure = unzip_this;

    success = 1;
}
void polymer::unzip_update(int s_index, int side){
    helix_list[s_index]->shorten(zip_unzip_structure);

    //update running centre
    int beta{ helix_list[s_index]->get_beta() };
    int kappa{ abs((1 - beta) - (1 - side)) };
    std::vector<double> rc{ helix_list[s_index]->get_rc(kappa) }, v{helix_list[s_index]->get_v()};
    double gamma{ pow(-1,kappa) };
    double pitch{ 0.28 };
    std::vector<double> temp{ multiplication_by_scalar(gamma*pitch,v) };
    rc = vector_addition(rc, temp);
    helix_list[s_index]->set_rc(kappa, rc);

    //change the structure id of the ones that got unzipped
    if (side == 0) {
        chain[zip_unzip_structure[0]]->change_structure(0);
        chain[zip_unzip_structure[3]]->change_structure(0);
    }
    else if (side == 1) {
        chain[zip_unzip_structure[1]]->change_structure(0);
        chain[zip_unzip_structure[2]]->change_structure(0);

    }


    update_positions();


}


std::vector<std::vector<double>> polymer::add_bp_to_helix(helix_struct* double_helix, int side) {
    int beta{ double_helix->get_beta() };
    int sigma{ side };
    std::vector<double> running_centre(3);
    sigma == 0 ? running_centre = double_helix->get_rc(beta) : running_centre = double_helix->get_rc(1-beta);

    std::vector<int> small_helix{ double_helix->get_monomers() };
    std::vector<double> r_x(3), r_y(3);

    int o_x, o_y;
    std::vector<int> extension{ small_helix };

    if (sigma == 0) {
        o_x = small_helix[0];
        o_y = small_helix[3];
        extension[0] = o_x - 1;
        extension[3] = o_y + 1;
    }
    else if (sigma == 1) {
        o_x = small_helix[1];
        o_y = small_helix[2];
        extension[1] = o_x + 1;
        extension[2] = o_y - 1;

    }
    r_x = chain[o_x]->get_position();
    r_y = chain[o_y]->get_position();

    std::vector<double> u(3), w(3);
    u = vector_subtraction(r_x,running_centre);
    //normalize(u);
    w = vector_subtraction(r_y, running_centre);
    //normalize(w);

    double phi = 33 * (atan(1) * 4) / 180;
    double radius{ 0.25 };
    double pitch{ 0.28 };

    std::vector<double> temp(3);
    std::vector<double> v{ double_helix->get_v() };

    if (abs((1-beta)-(1-sigma))== 0) {
        pitch = -pitch;
        phi = -phi;
    }

    temp = multiplication_by_scalar(pitch, v);
    running_centre = vector_addition(running_centre, temp);
    u = rotate_helix(u, v, phi);
    w = rotate_helix(w, v, phi);

    std::vector<double> n_x(3), n_y(3);
    n_x = vector_addition(running_centre, u);
    n_y = vector_addition(running_centre, w);


    return { n_x,n_y };

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// our double helix structures are represented by [a b c d] where a b c and d are not consecutive to each other usually, they
//define limits. sometimes its useful to "unpack" that representation: [a,a+1,...,b, c,c+1,...,d] and vice versa.
void polymer::unpack_region(std::vector<int>& reaction_region)
{
    if (reaction_region.size() != 4) { 
        std::cout << "tried to unpack region which was already unpacked" << std::endl;
        return; }

    int size{ abs(reaction_region[1] - reaction_region[0]) + 1 };
    std::vector<int> temp(2 * size);

    for (int i{ 0 }; i < size; i++) {
        temp[i] = reaction_region[0] + i;
    }
    for (int i{ 0 }; i < size; i++) {
        temp[i + size] = reaction_region[2] + i;
    }
    reaction_region = temp;

}
// our double helix structures are represented by [a b c d] where a b c and d are not consecutive to each other usually, they
//define limits. sometimes its useful to "unpack" that representation: [a,a+1,...,b, c,c+1,...,d] and vice versa.
void polymer::pack_region(std::vector<int>& reaction_region) {
    if (reaction_region.size() == 4) { 
        std::cout << "tried to pack region which was already packed" << std::endl;
        return; }

    int size{ static_cast<int>(reaction_region.size()) }, mid3p{ size / 2 - 1 }, mid5p{ mid3p + 1 };

    std::vector<int> temp(4);
    temp[0] = reaction_region[0];
    temp[1] = reaction_region[mid3p], temp[2] = reaction_region[mid5p];
    temp[3] = reaction_region[size - 1];

    reaction_region = temp;
}

void polymer::neighbouring_linkers()
{
    std::vector<std::vector<double>> np(N);
    for (auto i{ 0 }; i < N; i++) {
        np[i] = chain[i]->get_position();
    }
    new_config_positions = np;

    //std::vector<double> ws{ old_weights->get_weights() }, us{ old_weights->get_energies() };
    //new_weights->set_energies(us), new_weights->set_weights(ws);

    //std::cout << "neighbouring linkers called" << std::endl;
    std::vector<int> s(4);
    int limits;
    if (helix_list.size() > 0) {
        std::vector<int> linked_monomers;
        get_linked_monomers(linked_monomers);
        //print_1d_int_vec(linked_monomers);
        std::sort(linked_monomers.begin(), linked_monomers.end());

        for (int i{ 0 }; i < chain.size(); i++) {

            if (i == 0) {
                
                chain[i]->left_right_link[1] = linked_monomers.front();
            }
            else if (i == chain.size() - 1) {
                chain[i]->left_right_link[0] = linked_monomers.back();

            }
            else {
                //search left
                for (int l{ 1 }; l < i - 1; l++) {
                    if (chain[i - l]->get_structure() != 0) {// when you get to 
                        chain[i]->left_right_link[0] = i - l;
                        break;
                    }
                }
                //search right
                for (int r{ 1 }; r < chain.size() - i - 1; r++) {
                    if (chain[i + r]->get_structure() != 0) {
                        chain[i]->left_right_link[1] = i + r;
                        break;
                    }
                }
            }
        }

    }
    else {
        for (int i{ 0 }; i < chain.size(); i++) {
            chain[i]->left_right_link[0] = 0;
            chain[i]->left_right_link[1] = chain.size() - 1;
        }
    }

}

void polymer::update_large_struct_list() {
    zipped_structures.clear();

    for (int i{ 0 }; i < helix_list.size(); i++) {
        if (helix_list[i]->extended()) {
            zipped_structures.push_back(i);
        }
    }

}

void polymer::update_extensible_structures() {

    std::vector<int> temp;
    extendable_structures.clear();
    for (int s{ 0 }; s < helix_list.size(); s++) {
        if (!helix_list[s]->fully_extended()) {
            temp = helix_list[s]->get_monomers();

            int init_size{ abs(temp[1] - temp[0]) + 1 };
            structure_extension(temp, init_size);
            if (temp[1] - temp[0] + 1 > init_size) {
                print_1d_int_vec(temp);

                //pack_region(temp);
                //helix_list[s]->set_extendable(true);
                helix_list[s]->set_extension(temp);
                extendable_structures.push_back(s);
            }

        }
    }
}

void polymer::update_excluded_volume(std::vector<std::vector<int>>& growth_limits, std::vector<int> helix)// use this function once we've calculated the growth limits
//for a certain move. we update the excluded volume by adding all the monomers that will not be regrown to the excluded volume.
{
    excluded_volume.clear();
    std::vector<int> indices(chain.size());

    // we make a vector of indices which counts upwards from 0 to chain size. we will then delete the elements corresponding to monomers which will be regrown
    for (int i{ 0 }; i < chain.size(); i++) {
        indices[i] = i;
    }
    // go through all the limits and remove the monomers which will be regrown from the indices vector
    for (auto i : growth_limits) {
        std::sort(i.begin(), i.end());
        if (i[0] == 0) {
            //indices.erase(indices.begin(), indices.begin() + i[1]);
            std::fill(indices.begin()+i[0], indices.begin() + i[1], -1);
        }
        else if (i[1] == chain.size() - 1) {
            std::fill(indices.begin() + i[0]+1, indices.begin() + i[1]+1, -1);

        }
        else {
            //indices.erase(indices.begin() + i[0] + 1, indices.begin() + i[1]);
            std::fill(indices.begin() + i[0] + 1, indices.begin() + i[1], -1);

        }
    }
    // if the regrowth involves a helix (ie link) then we remove the elements of the monomers in the helix
    if (!helix.empty()) {
        unpack_region(helix);

        for (auto i : helix) {

            //indices.erase(indices.begin() + i);
            indices[i] = -1;
        }
        pack_region(helix);

    }

    // also, unless one of the growth limits is an endpoint (in which case it is not fixed) we should 
    // include the limit extrema of our growth in the excluded volume since they are fixed
    //for (auto i : growth_limits) {
    //    for(auto j:i){
    //        if(j != 0 && j != chain.size()-1){
    //            indices[j] == j;
    //        }
    //    }
    //}

    // finally we load the positions of the remaining monomers into the excluded volume 2d vector.
    std::vector<std::vector<double>> ev(indices.size() - 1);
    for (int j{ 0 }; j < indices.size() - 1; j++) {

        if (indices[j] != -1) {
            ev[j] = new_config_positions[indices[j]];

            //ev[j] = chain[indices[j]]->get_position();

        }
    }
    ev.erase(std::remove_if(ev.begin(), ev.end(),[](const std::vector<double>& innerVec) {return innerVec.empty();}),ev.end());
    excluded_volume = ev;
}

void polymer::helix_excluded_volume_interaction(std::vector<std::vector<double>>& helix_positions)
{
    double energy{ 0 };
    // we want to take into account the excluded volume interaction of the helix with the monomers that are not being regrown.
    if (excluded_volume.size() == 0) {
        helix_interaction_weight = 1;
    }
    else {
        for (int i{ 0 }; i < helix_positions.size(); i++) {
            energy = energy + u_r(helix_positions[i], excluded_volume);
        }

        helix_interaction_weight = exp(-energy);
    }

}


void polymer::grow_limits(std::vector<std::vector<int>>& limits, int alpha, bool forward_move)
{
    // forward_move  == true, we are growing limits for the new configuration
    // forward_move  == false, we are growing limits for the old configuration
    rosenbluth_growth* weight_object;
    std::vector<std::vector<double>>* positions;
    std::vector<std::vector<double>> position_subset;
    if (forward_move == true) {
        weight_object = new_weights;
        positions = &new_config_positions;

    }
    else {
        weight_object = old_weights;
        positions = &old_config_positions;
    }
    //for (auto i : limits) {
    //    print_1d_int_vec(i);
    //}
    //std::cout << "the length of excluded volume at start is " << excluded_volume->size()<<std::endl;
    for (auto i : weight_object->get_weights()) {
        if (i == 0.0) {
            std::cout << "problem" << std::endl;
        }
    }
    //bug catching : yamakawa growth, 1st monomer in the limit is not loaded in properly
    for (auto i : *positions) {
        if (i.size() == 0) {
            std::cout << "Problem" << std::endl;
        }
    }

    std::vector<double> W_vals, U_vals;// W weights, U energies
    int m0, kappa;
    std::vector<int> limit(2);
    int deletion_counter{ 0 };

    for (size_t i = 0; i < limits.size(); i++)
    {

        limit = limits[i];
        //std::cout << "limit" << std::endl;
        //print_1d_int_vec(limit);

        if (limit[0] == limit[1]) {
            std::cout << "growth limits were the same monomer" << std::endl;

            limits.erase(limits.begin() + i - deletion_counter);
            deletion_counter++;
            continue;
        }
        if (abs(limit[1] - limit[0]) <= 2) {
            limits.erase(limits.begin() + i - deletion_counter);
            deletion_counter++;
            continue;
        }
        //print_1d_int_vec(limit);
        W_vals.clear();//W weights
        U_vals.clear();// U energies
        std::vector<std::vector<double>> linker(abs(limit[0] - limit[1]) + 1), accepted_positions;
        kappa = static_cast<int>(pow(-1.0, alpha));
        m0 = 0;
        std::vector<std::vector<double>> fp;

        if (forward_move == false) {
            accepted_positions = get_subpositions(limit);
            if (limit[0] == alpha * (chain.size() - 1) && limit[1] == (1 - alpha) * (chain.size() - 1)) {
                accepted_positions.erase(accepted_positions.begin());

            }
            else {
                if (limit[0] == 0 || limit[0] == chain.size() - 1) {
                    accepted_positions.erase(accepted_positions.begin());
                }
                if (limit[1] == 0 || limit[1] == chain.size() - 1) {
                    accepted_positions.pop_back();
                }

            }
        }
        if (accepted_positions.size() >= linker.size()) {
            std::cout << "stop" << std::endl;
        }
        if (abs(static_cast<int>(linker.size()) - static_cast<int>(accepted_positions.size()))== 2) {
            std::cout << "stop" << std::endl;
        }

        //if one of the limits is an endpoint then we do random walk growth
        if (limit[1 - alpha] == chain.size() - 1 || limit[alpha] == 0) {
            if (limit[0] == chain.size() - 1 || limit[0] == 0) {
                kappa = -kappa;
                m0 = 1;
            }
            m0 = limit[m0];
            //std::cout << "size " << positions->operator[](m0).size() << std::endl;
            std::vector<double> o1{ positions->operator[](m0) };
            if (o1.size() == 0) {
                std::cout << "stop " << std::endl;
            }
            fp.push_back(o1);
            if (rosenbluth_switch == true) {
                if (forward_move == true) {
                    linker = grow_section(fp, linker.size() - 1, excluded_volume, U_vals, W_vals, false, forward_move);

                }
                else if(forward_move == false) {
                    linker = grow_section(fp, linker.size() - 1, excluded_volume, U_vals, W_vals, false, forward_move, accepted_positions);
                }
                for (auto i : linker) {
                    if (i.size() == 0) {
                        std::cout << "Problem" << std::endl;
                    }
                }

                //linker = rosenbluth_random_walk(o1, linker.size() - 1, *excluded_volume, U_vals, W_vals);

                // this is so that we don't include the fixed endpoint. 
                std::vector<int> temp_lim{ limit };
                if (limit[0] == 0) {
                    temp_lim[1] = limit[1] - 1;
                }
                if (limit[0] == chain.size() - 1) {
                    temp_lim[1] = limit[1] + 1;
                }
                if (limit[1] == 0) {
                    temp_lim[0] = limit[0] - 1;
                }
                if (limit[1] == chain.size() - 1) {
                    temp_lim[0] = limit[0] + 1;
                }

                weight_object->modify_weights(temp_lim, W_vals);
                weight_object->modify_energies(temp_lim, U_vals);


            }
            else if (rosenbluth_switch == false) {
                linker = random_walk(o1, linker.size() - 1, excluded_volume);
                W_vals.resize(linker.size());
                U_vals.resize(linker.size());

                std::fill(W_vals.begin(), W_vals.end(), 1.0);
                std::fill(U_vals.begin(), U_vals.end(), 0.0);
                weight_object->modify_weights(limit, W_vals);
                weight_object->modify_energies(limit, U_vals);

            }
            for (auto i : W_vals) {
                if (i == 0.0) {
                    std::cout << "Problem" << std::endl;
                }
            }

            int m;
            for (int i{ 0 }; i < linker.size(); i++) {
                m = m0 + kappa * i;
                positions->operator[](m) = linker[i];

            }

        }
        else {// if neither of the limits are endpoints then we do fixed endpoints growth
            //std::vector<double> p1{ chain[limit[0]]->get_position() }, p2{ chain[limit[1]]->get_position() };
            std::vector<double> p1{ positions->operator[](limit[0]) }, p2{ positions->operator[](limit[1]) };
            fp = { p1,p2 };
            if (rosenbluth_switch == true) {
                if (forward_move == true) {
                    linker = grow_section(fp, linker.size() - 1, excluded_volume, U_vals, W_vals, true, forward_move);

                }
                else if (forward_move == false) {
                    linker = grow_section(fp, linker.size() - 1, excluded_volume, U_vals, W_vals, true, forward_move, accepted_positions);
                }

                //linker = grow_section(fp, linker.size() - 1, *excluded_volume, U_vals, W_vals, true, forward_move);
                //linker = rosenbluth_grow_chain(p1, p2, linker.size() - 1, *excluded_volume, U_vals, W_vals);

                // this is so that we don't include the fixed endpoint.
                std::vector<int> temp_lim{ limit };
                if (alpha == 0) {
                    temp_lim[0] += 1;
                    temp_lim[1] -= 1;
                }
                else {
                    temp_lim[0] -= 1;
                    temp_lim[1] += 1;
                }
                weight_object->modify_weights(temp_lim, W_vals);
                weight_object->modify_energies(temp_lim, U_vals);




            }
            else if (rosenbluth_switch == false) {
                linker = grow_chain(p1, p2, linker.size() - 1, excluded_volume);

                // when we don't have rosenbluth sampling, instead of modifying other parts of code greatly, we will just set all the 
                // rosenbluth related factors to 1, that way they don't impact the result of the acceptance calculation (however with a 
                // cleaner solution we would save computation time).
                W_vals.resize(linker.size());
                U_vals.resize(linker.size());
                std::fill(W_vals.begin(), W_vals.end(), 1.0);
                std::fill(U_vals.begin(), U_vals.end(), 0.0);
                weight_object->modify_weights(limit, W_vals);
                weight_object->modify_energies(limit, U_vals);
            }
            for (int i{ 0 }; i < linker.size(); i++) {
                //chain[limit[0] + kappa * i]->change_position(linker[i]);
                positions->operator[](limit[0] + kappa * i) = linker[i];
            }

        }
        for (auto i : linker) {
            if (i.size() == 0) {
                std::cout << "Problem" << std::endl;
            }
        }

        //print_monomer_positions();
    }

    //std::cout << "the length of excluded volume at end is " << excluded_volume->size() << std::endl;
    //std::cout << "there are " << chain.size() << "monomers in total " << std::endl;

    for (auto i : weight_object->get_weights()) {
        if (i == 0.0) {
            std::cout << "problem" << std::endl;
        }
    }


    //bug catching : yamakawa growth, 1st monomer in the limit is not loaded in properly
    for (auto i : *positions) {
        if (i.size() == 0) {
            std::cout << "Problem" << std::endl;
        }
    }


    // sometimes we get limits which don't make much sense like [45,45] or limits which don't require new growth like [43,44]. 
    // instead of compensating for that in the growth_limits() functions e.g link_growth_limits() or in the move functions e.g link()
    // we do it here. it's a choice, not sure if the best one or not.
    if (forward_move == true) {
        new_growth_lims = limits;
        if (alpha == 1) {
            for (auto i : new_growth_lims) {
                std::sort(i.begin(), i.end());
            }
        }
    }
    else if (forward_move == false) {
        old_growth_lims = limits;
        if (alpha == 1) {
            for (auto i : old_growth_lims) {
                std::sort(i.begin(), i.end());
            }
        }

    }


}




void polymer::update_positions() {
    for (int i{ 0 }; i < chain.size(); i++) {
        if (!new_config_positions[i].empty()) {
            chain[i]->change_position(new_config_positions[i]);
        }
    }
    old_config_positions = new_config_positions;
    return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void polymer::swivel(bool success) {

    // sample a structure from the extended structures. p_select = 1 / extended_structures.size()
    int structure_index{ static_cast<int>(rand2(0,helix_list.size())) };
    helix_struct* double_helix{ helix_list[structure_index] };
    std::vector<double> u{ double_helix->get_u() }, v{ double_helix->get_v() };
    //std::vector<double> u{ helix_list[structure_index]->get_u() }, v{ helix_list[structure_index]->get_v() };

    // decide to sample a new u or v. p=1/2
    std::vector<double> u_n(3), v_n(3);
    double rn{ rand2(0,1) };

    double tau{0.5};
    //rn = 0.5;
    if (rn < 0.33333) { // v swivel

        // sample an angle
        double theta{ rand2(0,tau * atan(1) * 2) };
        int kappa{ static_cast<int>(rand2(0,2)) };
        theta = theta * pow(-1, kappa);

        if (reject_spin(double_helix, theta)) {
            //std::cout << "v swivel rejected" << std::endl;
            success = 0;
        }
        else {
            //std::cout << " v swivel success" << std::endl;

            std::vector<double> rotation(3);
            spin(double_helix, theta, rotation);

            // calculate acceptance rule. WIP.

            //UPDATES

            //update running centres
            std::vector<double> rc1{ double_helix->get_rc(0) }, rc2{ double_helix->get_rc(1) };
            std::vector<double> com{ centre_of_mass(double_helix->get_monomers()) };
            rc1 = rotate_by_angle_about_general_axis(rc1, rotation, theta, com);
            rc2 = rotate_by_angle_about_general_axis(rc2, rotation, theta, com);
            double_helix->set_rc(0, rc1);
            double_helix->set_rc(1, rc2);

            //std::cout << "running centres" << std::endl;
            //print_1d_doub_vec(rc1);
            //print_1d_doub_vec(rc2);

            //update v, helix vector
            std::vector<double> v_new{ vector_subtraction(rc2,rc1) };
            normalize(v_new);
            double_helix->change_v(v_new);

            //update positions
            update_positions();

            success = 1;
            return;
        }

    }
    else if (rn < 0.66666) { // u swivel

        // sample an angle
        double theta{ rand2(0,tau * atan(1) * 8)};
        int kappa{ static_cast<int>(rand2(0,2)) };
        theta = theta * pow(-1, kappa);

        if (reject_corkscrew(double_helix, theta)) {
            std::cout << "u swivel rejected" << std::endl;

            success = 0;
        }
        else {
            std::cout << " u swivel success" << std::endl;
            std::cout << "angle " << theta * 180 / (4 * atan(1)) << std::endl;

            corkscrew(double_helix, theta);

            // calculate acceptance rule. WIP.

            // since u is something we calculate not store, we don't need to update it. additionally, 
            // the running centres don't change since we just rotated around them. so just need to change 
            // the positions
            update_positions();

            success = 1;
        }
        

    }
    else { // origin translation
        double translation_distance{ 0.5 };
        std::vector<double> translation(3);
        sample_jump_direction(translation, translation_distance);

        if (reject_translate(double_helix, translation)) {
            std::cout << "o translation rejected" << std::endl;

            success = 0;
        }
        else {
            std::cout << "o translation success" << std::endl;

            translate(double_helix, translation);

            // calculate acceptance rule

            // if accepted, need to update positions, including the position of the running centre
            std::vector<double> rc1{ double_helix->get_rc(0) }, rc2{ double_helix->get_rc(1) };
            rc1 = vector_addition(rc1, translation);
            rc2 = vector_addition(rc2, translation);
            double_helix->set_rc(0, rc1);
            double_helix->set_rc(1, rc2);
            update_positions();
            success = 1;
        }
    }

}

bool polymer::reject_spin(helix_struct* double_helix, double angle) {
    std::vector<int> s{ double_helix->get_monomers() };
    std::vector<double> com{ centre_of_mass(s) };

    // calculate u
    int beta{ double_helix->get_beta() }, alpha{ double_helix->get_alpha() };
    std::vector<double> running_centre(3);
    beta == 0 ? running_centre = double_helix->get_rc(beta) : running_centre = double_helix->get_rc(1 - beta);
    int origin{ s[3 * alpha + (pow(-1,alpha) * beta)] };
    std::reverse(s.begin(), s.end());
    int origin_link{ s[3 * alpha + (pow(-1,alpha) * beta)] };
    std::reverse(s.begin(), s.end());
    std::vector<double> r_ox(3);
    r_ox = chain[origin]->get_position();
    std::vector<double> u(3);
    u = vector_subtraction(r_ox, running_centre);// our rotation axis

    // should make sure that the neighbouring linkers were already linked at this point

    // need to rotate the monomers defining a b c d of the structure and see if they maintain overstretching constraints
    if (alpha == 1) {
        std::sort(s.begin(), s.end(), std::greater<>());

    }
    int a{ s[0] }, b{ s[1] }, c{ s[2] }, d{ s[3] };
    if (alpha == 1) {
        std::sort(s.begin(), s.end());

    }

    int w, x, y, z;
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    std::vector<double> temp(3),temp2(3);

    // check w extreme end.
    if (w != alpha * (chain.size() - 1))
    {
        // rotate a 
        temp = chain[a]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, u, angle, com);
        
        temp2 = chain[w]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(w - a)) {
            return true;
        }
    }

    // z extreme end
    if (z != (1 - alpha) * (chain.size() - 1)) {

        //rotate d
        temp = chain[d]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, u, angle, com);

        temp2 = chain[z]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(d - z)) {
            return true;
        }
    }

    //check middle (enclosed) limits
    if (abs(b - c) > abs(b - x)) {

        //rotate b
        temp = chain[b]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, u, angle, com);

        temp2 = chain[x]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(b - x)) {
            return true;
        }
        
        temp = chain[c]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, u, angle, com);

        temp2 = chain[y]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(y - c)) {
            return true;
        }
        
    }
    return false;

}

bool polymer::reject_corkscrew(helix_struct* dh, double angle) {

    std::vector<int> s{ dh->get_monomers() };
    std::vector<double> v{ dh->get_v() };
    normalize(v);

    int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };
    std::vector<double> running_centre1(3), running_centre2(3);
    if (alpha == 1) {
        std::sort(s.begin(), s.end(), std::greater<>());

    }
    int a{s[0]},b{s[1]},c{s[2]},d{s[3]};
    if (alpha == 1) {
        std::sort(s.begin(), s.end());

    }

    int w, x, y, z;
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    if (beta == 0) {
        running_centre1 = dh->get_rc(beta), running_centre2 = dh->get_rc(1 - beta);

    }
    else {
        running_centre1 = dh->get_rc(1 - beta), running_centre2 = dh->get_rc(beta);

    }

    std::vector<double> temp(3), temp2(3);

    // check w extreme end.
    if (w != alpha * (chain.size() - 1))
    {
        // rotate a 
        temp = chain[a]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, v, angle, running_centre1);

        temp2 = chain[w]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(w - a)) {
            return true;
        }
    }

    // z extreme end
    if (z != (1 - alpha) * (chain.size() - 1)) {

        //rotate d
        temp = chain[d]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, v, angle, running_centre1);

        temp2 = chain[z]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(d - z)) {
            return true;
        }
    }

    //check middle (enclosed) limits
    if (abs(b - c) > abs(b - x)) {

        //rotate b
        temp = chain[b]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, v, angle, running_centre2);

        temp2 = chain[x]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(b - x)) {
            return true;
        }

        temp = chain[c]->get_position();
        temp = rotate_by_angle_about_general_axis(temp, v, angle, running_centre2);

        temp2 = chain[y]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(y - c)) {
            return true;
        }

    }
    return false;


}

void polymer::spin(helix_struct* dh, double angle, std::vector<double>& rotation_axis) {
    //helix_struct* dh{ sample_double_helix() };
    std::vector<int> s{ dh->get_monomers() };

    std::vector<double> com{ centre_of_mass(s) };
    
    int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };
    std::vector<double> v_p{ dh->get_v() };
    std::vector<double> rot_axis{ sample_unit_perp_vector(v_p) };
    rotation_axis = rot_axis;

    unpack_region(s);
    std::vector<std::vector<double>> rotated_positions(s.size());
    std::vector<double> temp(3);
    for (int i{ 0 }; i < s.size();i++) {
        temp = chain[s[i]]->get_position();
        rotated_positions[i] = rotate_by_angle_about_general_axis(temp, rot_axis, angle, com);
        new_config_positions[s[i]] = rotated_positions[i];
    }

    // should test that with our calculated v we do indeed generate the desired helix.
    //std::cout << "new v vector calculated by rotating running centres and taking their vector difference" << std::endl;
    std::vector<double> rc1{ dh->get_rc(0) }, rc2{ dh->get_rc(1) };
    rc1 = rotate_by_angle_about_general_axis(rc1, rot_axis, angle, com);
    rc2 = rotate_by_angle_about_general_axis(rc2, rot_axis, angle, com);
    pack_region(s);
    std::vector<double> temp1(1), temp2(3);
    if (beta == 0) {
        temp1 = new_config_positions[s[0]];
        temp2 = new_config_positions[s[3]];

        temp1 = new_config_positions[s[1]];
        temp2 = new_config_positions[s[2]];

    }
    else if (beta == 1) {
        temp1 = new_config_positions[s[0]];
        temp2 = new_config_positions[s[3]];

        temp1 = new_config_positions[s[1]];
        temp2 = new_config_positions[s[2]];

    }
    //std::cout << "v theory 0" << std::endl;
    std::vector<double> v_theory{ vector_subtraction(rc1,rc2) };
    normalize(v_theory);
    //print_1d_doub_vec(v_theory);


    std::vector<std::vector<int>> growth_limits;
    swivel_growth_limits(dh, growth_limits);
    update_excluded_volume(growth_limits);

    grow_limits(growth_limits, alpha, true);

}

std::vector<double> polymer::centre_of_mass(std::vector<int> s) {
    std::vector<double> r_com{chain[s[0]]->get_position()};
    unpack_region(s);
    std::vector<double> temp(3);
    for (int i{ 1 }; i < s.size(); i++) {
        temp = chain[s[i]]->get_position();
        r_com = vector_addition(r_com, temp);
    }
    int k{ static_cast<int>(s.size()) };
    std::for_each(r_com.begin(), r_com.end(), [k](double& c) { c /= k; });
    //std::cout << "centre of mass vector " << std::endl;
    //print_1d_doub_vec(r_com);
    return r_com;
}

void polymer::corkscrew_update(helix_struct* dh)
{
    update_positions();


}

void polymer::spin_update(helix_struct* dh, std::vector<double>& rot_axis, double angle)
{
    std::vector<double> rc1{ dh->get_rc(0) }, rc2{ dh->get_rc(1) };
    std::vector<double> com{ centre_of_mass(dh->get_monomers()) };
    rc1 = rotate_by_angle_about_general_axis(rc1, rot_axis, angle, com);
    rc2 = rotate_by_angle_about_general_axis(rc2, rot_axis, angle, com);
    dh->set_rc(0, rc1);
    dh->set_rc(1, rc2);

    //update v, helix vector
    std::vector<double> v_new{ vector_subtraction(rc2,rc1) };
    normalize(v_new);
    dh->change_v(v_new);

    //update positions
    update_positions();

}

void polymer::translate_update(helix_struct* dh, std::vector<double>& translation)
{
    // if accepted, need to update positions, including the position of the running centre
    std::vector<double> rc1{ dh->get_rc(0) }, rc2{ dh->get_rc(1) };
    rc1 = vector_addition(rc1, translation);
    rc2 = vector_addition(rc2, translation);
    dh->set_rc(0, rc1);
    dh->set_rc(1, rc2);
    update_positions();

}

double polymer::configuration_energy()
{
    std::vector<std::vector<double>> positions(chain.size());

    positions = generate_vector_of_monomer_positions();

    double cumulative_energy{ 0 };
    for (auto i : positions) {
        cumulative_energy += u_r(i, positions);
    }
    return cumulative_energy;
}

void polymer::corkscrew(helix_struct* dh, double angle) {

    std::vector<int> s{ dh->get_monomers() };
    std::vector<double> v{ dh->get_v() };
    normalize(v);

    int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };
    std::vector<double> running_centre(3);
    beta == 0 ? running_centre = dh->get_rc(beta) : running_centre = dh->get_rc(1 - beta);

    int n{ abs(s[0] - s[1]) + 1 };
    std::vector<std::vector<double>> rotated_positions(2*n);

    // starting pair to be rotated depends on beta
    int m1, m2;
    if (beta == 0)
    {
        m1 = s[0];
        m2 = s[3];
    }
    else {
        m1 = s[1];
        m2 = s[2];
    }
    int gamma{ static_cast<int>(pow(-1,beta)) };
    double pitch{ 0.28 };
    std::vector<double> temp{ multiplication_by_scalar(pitch, v) };

    std::vector<double> p1{ chain[m1]->get_position() }, p2{ chain[m2]->get_position() };

    p1 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);
    p2 = rotate_by_angle_about_general_axis(p2, v, angle, running_centre);

    new_config_positions[m1] = p1, new_config_positions[m2] = p2;
    //chain[m1]->change_position(p1), chain[m2]->change_position(p2);
    for (int i = 1; i < n; i++)
    {
        running_centre = vector_addition(running_centre, temp);

        m1 = m1 + gamma;
        m2 = m2 - gamma;

        p1 = chain[m1]->get_position();
        p2 = chain[m2]->get_position();

        p1 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);
        p2 = rotate_by_angle_about_general_axis(p2, v, angle, running_centre);

        new_config_positions[m1] = p1, new_config_positions[m2] = p2;
        //chain[m1]->change_position(p1), chain[m2]->change_position(p2);


    }

    std::vector<std::vector<int>> growth_limits;
    swivel_growth_limits(dh, growth_limits);
    update_excluded_volume(growth_limits);

    grow_limits(growth_limits, alpha, true);

}

bool polymer::reject_translate(helix_struct* dh, std::vector<double>& translation) {
    std::vector<int> s{ dh->get_monomers() };
    int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };

    if (alpha == 1) {
        std::sort(s.begin(), s.end(), std::greater<>());

    }
    int a{ s[0] }, b{ s[1] }, c{ s[2] }, d{ s[3] };
    if (alpha == 1) {
        std::sort(s.begin(), s.end());

    }

    int w, x, y, z;
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    std::vector<double> temp(3), temp2(3);

    // check w extreme end.
    if (w != alpha * (chain.size() - 1))
    {
        // rotate a 
        temp = chain[a]->get_position();
        temp = vector_addition(temp,translation);

        temp2 = chain[w]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(w - a)) {
            return true;
        }
    }

    // z extreme end
    if (z != (1 - alpha) * (chain.size() - 1)) {

        //rotate d
        temp = chain[d]->get_position();
        temp = vector_addition(temp, translation);

        temp2 = chain[z]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(d - z)) {
            return true;
        }
    }

    //check middle (enclosed) limits
    if (abs(b - c) > abs(b - x)) {

        //rotate b
        temp = chain[b]->get_position();
        temp = vector_addition(temp, translation);

        temp2 = chain[x]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(b - x)) {
            return true;
        }

        temp = chain[c]->get_position();
        temp = vector_addition(temp, translation);

        temp2 = chain[y]->get_position();
        if (dist_2_points3d(temp, temp2) > abs(y - c)) {
            return true;
        }

    }
    return false;

    
}

void polymer::translate(helix_struct* dh, std::vector<double> t) {
    std::vector<int> s{ dh->get_monomers() };
    int alpha{ dh->get_alpha() };
    unpack_region(s);

    std::vector<double> temp(3);
    for (auto i : s) {
        temp = chain[i]->get_position();
        new_config_positions[i] = vector_addition(temp, t);
    }
    std::vector<std::vector<int>> growth_limits;
    swivel_growth_limits(dh, growth_limits);
    update_excluded_volume(growth_limits);
    grow_limits(growth_limits, alpha, true);


}

void polymer::swivel_growth_limits(helix_struct* dh, std::vector<std::vector<int>>& limits) {
    std::vector<int> s{ dh->get_monomers() };
    int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };

    if (alpha == 1) {
        std::sort(s.begin(), s.end(), std::greater<>());

    }
    int a{ s[0] }, b{ s[1] }, c{ s[2] }, d{ s[3] };
    if (alpha == 1) {
        std::sort(s.begin(), s.end());

    }

    int w, x, y, z;
    w = chain[a]->left_right_link[alpha];
    x = chain[b]->left_right_link[1 - alpha];
    y = chain[c]->left_right_link[alpha];
    z = chain[d]->left_right_link[1 - alpha];

    std::vector<std::vector<int>> growth_limits;

    growth_limits.push_back({ w,a });
    if (abs(b - c) > abs(b - x)) {
        growth_limits.push_back({ b,x });
        growth_limits.push_back({ y,c });
    }
    else {
        growth_limits.push_back({ b,c });
    }
    growth_limits.push_back({ d,z });

    limits = growth_limits;

}