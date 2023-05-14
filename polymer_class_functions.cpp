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
    accepted_weights = new rosenbluth_growth(limit, weights, energies, 1);
    proposed_weights = accepted_weights;

    regrown_positions = positions;
}

polymer::polymer(std::vector<double> beginning_position, std::vector<double> end_position) {
    N = 2;
    chain.push_back(new monomer(beginning_position, 0, "U"));
    chain.push_back(new monomer(end_position, 1, "G"));
}

polymer::polymer(std::vector<std::vector<double>> generated_positions)
{
    std::vector<monomer*> M(generated_positions.size());
    for (int i{ 0 }; i < generated_positions.size(); i++) {
        M[i] = new monomer(generated_positions[i], i, random_base());
        //M.push_back(new monomer(generated_positions[i], i, random_base()));
        M[i]->left_right_link[1] = generated_positions.size()-1;
        //M[i]->get_position();
        //std::cout << M[i]->get_base() << " ";
    }
    //std::cout << std::endl;
    chain = M;
    search_results = structure_search();
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
int polymer::chain_length() {
    //std::cout << chain.size() << "\n";
    return chain.size();
}

void polymer::get_monomer_position(int id) {
    chain[id]->get_position();
}

void polymer::add_monomer(std::vector<double> position, std::string base_type) {
    auto iterator = chain.end() - 1;
    chain[chain.size() - 1]->id++;
    chain.insert(iterator, new monomer(position, chain.size() - 1, base_type));

    N++;
}

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
std::vector<std::string> particle_types{ "As","Gd","Hf","Mg","K", "P"};
void polymer::output_for_ovito6(std::string filename) {

    std::vector<int> linked_monomers;
    get_linked_monomers(linked_monomers);
    std::vector<std::vector<int>> structures;
    std::vector<int> s;
    for (auto i : helix_list) {
        s = i->get_monomers();
        structures.push_back(s);
    }
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

monomer* polymer::operator[](int i) {
    return chain[i];
}
//**************************************************************************************//
//**************************************************************************************//
//**************************************************************************************//


//**************************************************************************************//
//**************************************************************************************//
//**************************************************************************************//
bool polymer::compatible_bases(monomer* monomer_i, monomer* monomer_j) {
    if (monomer_i->get_base() == "U" && monomer_j->get_base() == "A" || monomer_i->get_base() == "A" && monomer_j->get_base() == "U") { return true; }
    else if (monomer_i->get_base() == "G" && monomer_j->get_base() == "C" || monomer_i->get_base() == "C" && monomer_j->get_base() == "G") { return true; }
    else { return false; }
}

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

std::vector<int> polymer::sample_structure(std::vector<std::vector<std::vector<int>>> potential_structures)
{
    int r{ 1 };
    //int r{ static_cast<int>(rand2(0,potential_structures.size())) };
    int k{ static_cast<int>(rand2(0,potential_structures[r].size())) };

    //std::cout << "sampled structure is " << std::endl;
    //for (auto i : potential_structures[r][k]) {
    //    std::cout << i << " ";
    //}std::cout << std::endl;

    return potential_structures[r][k];
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
    if (search_results.size() == 0) {
        link = {};
        return;
    }

    //// we're only considering regions of length 3 for the moment. the default argument of search results function is for r=3.
    int region_index0{ static_cast<int>(rand2(0,search_results.size())) };
    ss_index = region_index0;
    link = search_results[region_index0];

    // sample helix origin. 
    alpha = static_cast<int>(rand2(0, 2));
    beta = static_cast<int>(rand2(0, 2));

    std::cout << "sample link region done" << std::endl;

}
helix_struct* polymer::sample_double_helix() {
    //pick structure to unlink and an origin within the structure to start the 'degrowth'
    int s_i{ static_cast<int>(rand2(0,helix_list.size())) };
    return helix_list[s_i];

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


void polymer::unlink() {
    if (helix_list.size() == 0) {
        std::cout << "Cannot perform an unlink move" << std::endl;
    }
    else {
        neighbouring_linkers();

        //pick structure to unlink and an origin within the structure to start the 'degrowth'
        int s_i{ static_cast<int>(rand2(0,helix_list.size())) };
        std::vector<int> destructure{ helix_list[s_i]->get_monomers() };

        int alpha{ static_cast<int>(rand2(0, 2)) };//alpha == 0 is 3' to 5' growth. alpha == 1 is 5' to 3' growth
        int beta{ static_cast<int>(rand2(0, 2)) };
        std::cout << "alpha = " << alpha << " beta = " << beta << std::endl;

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

        std::cout << "abcd " << a << " " << b << " " << c << " " << d << std::endl;
        std::cout << "wxyz " << w << " " << x << " " << y << " " << z << std::endl;

        int region_size{ abs(b - a) + 1 };
        bool middle_struct{ true };
        if (abs(b - d) < abs(b - x)) {
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
            degrowth_limits.push_back({ w,b });
            degrowth_limits.push_back({ y,z });
        }

        std::cout << "degrowth limits" << std::endl;
        for (auto i : degrowth_limits) {
            std::cout << i[0] << " " << i[1] << std::endl;
        }


        //// perform degrowth
        // the growth bit. 
        update_excluded_volume(degrowth_limits);
        grow_limits(degrowth_limits, alpha);

        delete helix_list[s_i];
        helix_list.erase(helix_list.begin() + s_i);
        search_results.push_back(destructure); // since we are unlinking this region, it can now return to the pool of structures which can be made via a link move.

        //update neighbours
        neighbouring_linkers();

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

    // update the unlinked monomers. they are not part of a structure anymore.
    unpack_region(s);
    for (auto i : s) {
        chain[i]->change_structure(0);
    }

    //delete corresponding helix object and delete element from array.
    delete helix_list[s_i];
    helix_list.erase(helix_list.begin() + s_i);

    // update positions and neighbours
    update_positions();
    neighbouring_linkers();

    //accepted_weights = proposed_weights;

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
    neighbouring_linkers();//update neighbours
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
    std::vector<std::vector<int>> free_growth_limits;
    //std::cout << "abcd " << a << " " << b << " " << c << " " << d << std::endl;
    //std::cout << "wxyz " << w << " " << x << " " << y << " " << z << std::endl;

    if (beta == 1) {
        free_growth_limits.push_back({ w,a });
    }
    else if (beta == 0 && y != b) {
        free_growth_limits.push_back({ b,x });
    }
    free_growth_limits.push_back({ y,c });
    free_growth_limits.push_back({ d,z });



    //start the growth. 

    //first update the excluded volume
    //std::vector<std::vector<double>> excluded_volume;
    update_excluded_volume(free_growth_limits, reaction_region);

    // grow helix
    std::vector<std::vector<double>> s_x, s_y;
    std::vector<double> u, v;
    sample_helix_vectors(u, v);
    helix_struct* dh = new helix_struct(reaction_region,v,alpha,beta);
    print_1d_int_vec(reaction_region);

    std::vector<double> o_h{ chain[origin]->get_position() };
    std::vector<std::vector<double>> rcentres(2);
    generate(region_size, v, u, o_h, s_x, s_y,rcentres);
    dh->set_running_x_strand_centre(rcentres[0]);
    dh->set_running_y_strand_centre(rcentres[1]);

    int gamma{ static_cast<int>(pow(-1.0,alpha + beta)) };
    for (int i{ 0 }; i < region_size; i++) {
        regrown_positions[origin + gamma * i] = s_x[i];
        regrown_positions[origin_y + gamma * i] = s_y[i];

        //chain[origin + gamma * i]->change_position(s_x[i]);
        //chain[origin_y + gamma * i]->change_position(s_y[i]);
        excluded_volume.push_back(s_x[i]);
        excluded_volume.push_back(s_y[i]);
    }
    //print_monomer_positions();

    grow_limits(free_growth_limits, alpha);


    // update the linked neighbours for each monomer
    helix_list.push_back(dh);
    unpack_region(reaction_region);
    for (auto i : reaction_region) {
        chain[i]->change_structure(helix_list.size() + 1);
    }
    search_results.erase(search_results.begin() + ss_index);//this line ensures that we remove the linked structure from the pool of

    //search_results.erase(search_results.begin(), search_results.begin() + ss_index);//this line ensures that we remove the linked structure from the pool of


    ////trying to catch a bug where there's a huge separation
    //for (auto i : chain) {
    //    int j = i->get_id();
    //    if (j + 1 > chain.size() - 1) {
    //        break;
    //    }
    //    std::vector<double> temp1{ chain[j + 1]->get_position() }, temp2{ i->get_position() };
    //    if (dist_2_points3d(temp1, temp2) > 5) {
    //        std::cout << "CAUGHT ONE" << std::endl;
    //        std::cout << "monomer " << i->get_id() << " and " << j + 1 << std::endl;
    //    }
    //}

    neighbouring_linkers();//update neighbours


}


bool polymer::reject_link1(std::vector<int>& link_region, int alpha, int beta)
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

double polymer::acceptance_link0(bool l_or_u) {
    double W_b, W_u, Z_ni_A{ 0.27 }, Z_ni_B{ 0.27 }, Z_ni_AB{ 1 }, G0{ 0.5 }, rho0{ 1 }, U_b, U_u;
    double beta{ 1 };
    if (helix_list.size() == 0) {
        std::cout << "no link" << std::endl;
        return NULL;
    }
    else {
        helix_struct* link{ helix_list[0] };
        std::vector<int> s{ link->get_monomers() };
        int hairpin_segments{ abs(s[1] - s[2]) };
        double delta{ helix_separation() };
        double hairpin_gen_prob{ ideal_chain_pdf(delta,hairpin_segments) };

        W_u = hairpin_section->get_rosen_weight();
        U_u = hairpin_section->get_energy();

        W_b = free_section->subsection_weight({ s[0], s[3] });
        U_b = free_section->subsection_energy({ s[0], s[3] });

        double z = exp(beta * (U_u - U_b)) * (W_b / W_u) * (exp(-beta*G0) / rho0) * (Z_ni_A * Z_ni_B / Z_ni_AB) * hairpin_gen_prob;
        if (l_or_u == 1) {
            if (z < 1) {
                return z;
            }
            else { return 1.0; }

        }
        else if (l_or_u == 0) {
            if (1 / z < 1) {
                return 1 / z;
            }
            else {
                return 1.0;
            }
        }
    }

}

double polymer::regeneration_probabilities() {
    
    std::vector<double> old_probs(regrowth_lims.size()), new_probs(regrowth_lims.size());
    for (int i{ 0 }; i < regrowth_lims.size();i++) {
        old_probs[i] = accepted_weights->subsection_probability(regrowth_lims[i]);
        new_probs[i] = proposed_weights->subsection_probability(regrowth_lims[i]);
    }
    double old_config_prob{product_of_elements(old_probs)}, new_config_prob{product_of_elements(new_probs)};

    // are our limits ascending?
    // if we have two sections which are between two fixed endpoints then we need to use the probability distribution (yamakawa).
    std::vector<double> extra_prob(3);
    std::vector<int> limit(2);

    for (int i{ 0 }; i < regrowth_lims.size(); i++) {
        limit = regrowth_lims[i];
        if (limit[0] != 0 && limit[1] != chain.size() - 1) {
            std::vector<double> temp1(3),temp2(3);
            temp1 = chain[limit[0]]->get_position();
            temp2 = chain[limit[1]]->get_position();

            temp1 = vector_subtraction(temp1,temp2);
            double temp{ vector_modulus(temp1) };
            extra_prob[i] = ideal_chain_pdf(temp, abs(limit[1] - limit[0]));
        }
        else {
            extra_prob[i] = 1;
        }

    }
    return 0; 

}
void polymer::link_update(int ss_index, int a, int b, std::vector<double>&v_vec, std::vector<std::vector<double>> &rcentres) {
    std::vector<int> s{ search_results[ss_index] };

    helix_list.push_back(new helix_struct(s, v_vec, rcentres, a, b));

    unpack_region(s);
    for (auto i : s) {
        chain[i]->change_structure(helix_list.size() + 1);
    }
    search_results.erase(search_results.begin() + ss_index);//this line ensures that we remove the linked structure from the pool of

    neighbouring_linkers();//update neighbours

    //accepted_weights = proposed_weights;
    update_positions();

}

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

void polymer::sample_zip_move() {

}
void polymer::zip(bool& success) {

    neighbouring_linkers();//update neighbours

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

    sample_zip_region(s_z, reaction_region,side);
    print_1d_int_vec(s_z);

    if (s_z.size() == 0) {
        success = 0;
        return;
    }
    if (s_z[0]<0 || s_z[3]>chain.size() - 1) {
        std::cout << "invalid zip region" << std::endl;
        success = 0;
        return;
    }
    if (zip_structure_overlap(s_z, side)) {
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
    std::vector<double> temp(3);

    //////////////// what follows should be in its own function if the zip gets accepted ///////////////////

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

    temp = multiplication_by_scalar(pitch, v);
    running_centre = vector_addition(running_centre, temp);

    side == 0 ? extend_helix->set_rc(beta, running_centre) : extend_helix->set_rc(1 - beta, running_centre);

    //********************************************************

    // change the positions of the zipped monomers
    int kappa{ static_cast<int>(pow(-1,side)) };
    regrown_positions[o_x - kappa] = new_pos[0];
    regrown_positions[o_y + kappa] = new_pos[1];
    //chain[o_x - kappa]->change_position(new_pos[0]);
    excluded_volume.push_back(new_pos[0]);
    //chain[o_y + kappa]->change_position(new_pos[1]);
    excluded_volume.push_back(new_pos[1]);

    extend_helix->extend(s_z);

    //std::vector<std::vector<double>> excluded_volume;// will need to make a function which correctly gets excldued volume at the start
    //// as a function of growth limits.

    grow_limits(growth_limits, 0);

    //update the structure that we zipped. should also remove the one that's already been zipped.
    //helix_list[structure_index]->extend(s_z);
    unpack_region(s_z);
    for (auto i : s_z) {
        chain[i]->change_structure(structure_index + 1);
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
    std::cout << "sigma " << side << std::endl;

    std::cout << "ZIP MOVE SUCCESS" << std::endl;
    success = true;
    ////print_monomer_positions();
    //std::cout << "zip move success, region enlarged" << std::endl;
    //neighbouring_linkers();//update neighbours

}

void polymer::unzip(bool &success) {
    neighbouring_linkers();//update neighbours

    update_large_struct_list();
    if (zipped_structures.size() == 0) {
        std::cout << "unzip move not possible, no regions to unzip" << std::endl;
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
    std::cout << "unzip went through: region" << std::endl;
    print_1d_int_vec(unzip_this);
    std::vector<std::vector<int>> limits;
    unzip_growth_linits(zipped_helix, side, limits);

    update_excluded_volume(limits);
    grow_limits(limits, 0);

    //this section also needs to be changed once the acceptance probability is properly implemented
    // ie it needs to exist as an independent function which is called iff the move is not rejected 
    // due to the acceptance probability.
    helix_list[structure_index]->shorten(unzip_this);
    //change the structure id of the ones that got unzipped
    if (side == 0) {
        chain[zipped_helix[0]]->change_structure(0);
        chain[zipped_helix[3]]->change_structure(0);
    }
    else if (side == 1) {
        chain[zipped_helix[1]]->change_structure(0);
        chain[zipped_helix[2]]->change_structure(0);

    }
    //unpack_region(s_z);
    //structure_list[structure_index] = s_z;

    neighbouring_linkers();//update neighbours
    //output_for_ovito4("s3");
    //exit(0);
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
    std::vector<int> s(4);
    int limits;
    if (helix_list.size() > 0) {
        std::vector<int> linked_monomers;
        for (auto i : helix_list) {
            s = i->get_monomers();
            limits = s[1] - s[0];
            for (int j{ 0 }; j < limits; j++) {
                linked_monomers.push_back(s[0] + j);
                linked_monomers.push_back(s[2] + j);
            }
        }
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

    // finally we load the positions of the remaining monomers into the excluded volume 2d vector.
    std::vector<std::vector<double>> ev(indices.size() - 1);
    for (int j{ 0 }; j < indices.size() - 1; j++) {

        if (indices[j] != -1) {
            ev[j] = chain[indices[j]]->get_position();

        }
    }
    ev.erase(std::remove_if(ev.begin(), ev.end(),[](const std::vector<double>& innerVec) {return innerVec.empty();}),ev.end());
    excluded_volume = ev;
}

void polymer::grow_limits(std::vector<std::vector<int>>& limits, int alpha) {
   
    if (alpha == 1) {
        for (auto i : regrowth_lims) {
            std::sort(i.begin(), i.end());
        }
    }
    std::vector<double> W_vals, U_vals;
    int m0, kappa;
    std::vector<int> limit(2);
    for (size_t i = 0; i < limits.size(); i++)
    {
        limit = limits[i];

        if (limit[0] == limit[1]) {
            limits.erase(limits.begin() + i);
            continue;
        }
        //print_1d_int_vec(limit);
        W_vals.clear();
        U_vals.clear();
        std::vector<std::vector<double>> linker(abs(limit[0] - limit[1]) + 1);
        kappa = static_cast<int>(pow(-1.0, alpha));
        m0 = 0;

        //if one of the limits is an endpoint then we do random walk growth
        if (limit[1 - alpha] == chain.size() - 1 || limit[alpha] == 0) {
            if (limit[0] == chain.size() - 1 || limit[0] == 0) {
                kappa = -kappa;
                m0 = 1;
            }
            m0 = limit[m0];
            //linker = random_walk(chain[m0]->get_position(), linker.size() - 1, excluded_volume);

            //std::vector<double> o1{ chain[m0]->get_position() };
            std::vector<double> o1{ regrown_positions[m0] };

            linker = rosenbluth_random_walk(o1, linker.size() - 1, excluded_volume, U_vals, W_vals);
            //proposed_weights->modify_weights(limit, W_vals);
            //proposed_weights->modify_energies(limit, U_vals);


            int m;
            for (int i{ 0 }; i < linker.size(); i++) {
                m = m0 + kappa * i;
                regrown_positions[m] = linker[i];
                //chain[m]->change_position(linker[i]);
                
            }

        }
        else {// if neither of the limits are endpoints then we do fixed endpoints growth
            //std::vector<double> p1{ chain[limit[0]]->get_position() }, p2{ chain[limit[1]]->get_position() };
            std::vector<double> p1{ regrown_positions[limit[0]] }, p2{ regrown_positions[limit[1]] };

            //linker = grow_chain(chain[limit[0]]->get_position(), chain[limit[1]]->get_position(), linker.size() - 1, excluded_volume);
            linker = rosenbluth_grow_chain(p1, p2, linker.size() - 1, excluded_volume, U_vals,W_vals);
            //proposed_weights->modify_weights(limit, W_vals);
            //proposed_weights->modify_energies(limit, U_vals);

            for (int i{ 0 }; i < linker.size(); i++) {
                //chain[limit[0] + kappa * i]->change_position(linker[i]);
                regrown_positions[limit[0] + kappa * i] = linker[i];
            }

        }
        //print_monomer_positions();
    }
    regrowth_lims = limits;

}

void polymer::update_positions() {
    for (int i{ 0 }; i < chain.size(); i++) {
        if (!regrown_positions[i].empty()) {
            chain[i]->change_position(regrown_positions[i]);
        }
    }
    //std::vector<std::vector<double>> pos(N);
    //regrown_positions = pos;
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
    rn = 0.7;
    // these functions should check physical constraints, calculate acceptance values which will later be 
    // changed if the move is accepted.
    if (rn < 0.33333) { // v swivel

        // sample an angle
        double theta{ rand2(0,tau * atan(1) * 2) };
        int kappa{ static_cast<int>(rand2(0,2)) };
        theta = theta * pow(-1, kappa);

        if (reject_v_swivel(double_helix, theta)) {
            std::cout << "v swivel rejected" << std::endl;
            success = 0;
        }
        else {
            std::cout << " v swivel success" << std::endl;
            std::cout << "angle " << theta * 180 / (4 * atan(1)) << std::endl;
            v_rotate_double_helix(double_helix, theta);

            // calculate acceptance rule. WIP.

            // if accepted, this is the update. we have to change one of the running centres and v in our 
            // helix struct object . and we have to change the position.

            //following block is all to update the running centre and v.
            int beta{ double_helix->get_beta() };
            std::vector<double> running_centre(3);
            beta == 0 ? running_centre = double_helix->get_rc(beta) : running_centre = double_helix->get_rc(1 - beta);
            std::vector<double> old_v{ double_helix->get_v() };
            std::vector<double> new_v{ rotate_helix(old_v, u, theta) };
            double_helix->change_v(new_v);
            std::vector<int> s{ double_helix->get_monomers() };
            int n{ abs(s[1] - s[0])};
            double pitch{ 0.28 };
            std::vector<double> temp(3);
            temp = multiplication_by_scalar(n * pitch, new_v);
            running_centre = vector_addition(running_centre, temp);
            beta == 0 ? double_helix->set_rc(1 - beta, running_centre) : double_helix->set_rc(beta, running_centre);

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

        if (reject_u_swivel(double_helix, theta)) {
            std::cout << "u swivel rejected" << std::endl;

            success = 0;
        }
        else {
            std::cout << " u swivel success" << std::endl;
            std::cout << "angle " << theta * 180 / (4 * atan(1)) << std::endl;

            u_rotate_double_helix(double_helix, theta);

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

        if (reject_o_translation(double_helix, translation)) {
            std::cout << "o translation rejected" << std::endl;

            success = 0;
        }
        else {
            std::cout << "o translation success" << std::endl;

            origin_translation(double_helix, translation);

            // calculate acceptance rule

            // if accepted, need to update positions, nothing else.
            update_positions();
            success = 1;
        }
    }

}

bool polymer::reject_v_swivel(helix_struct* double_helix, double angle) {
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

bool polymer::reject_u_swivel(helix_struct* dh, double angle) {

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

void polymer::v_rotate_double_helix(helix_struct* dh, double angle) {
    //helix_struct* dh{ sample_double_helix() };
    std::vector<int> s{ dh->get_monomers() };

    std::vector<double> com{ centre_of_mass(s) };
    

    int beta{ dh->get_beta() }, alpha{ dh->get_alpha() };

    std::vector<double> running_centre(3);

    beta == 0 ? running_centre = dh->get_rc(beta) : running_centre = dh->get_rc(1 - beta);
    int origin{ s[3 * alpha + (pow(-1,alpha) * beta)] };
    std::reverse(s.begin(), s.end());

    int origin_link{ s[3 * alpha + (pow(-1,alpha) * beta)] };
    std::reverse(s.begin(), s.end());

    std::vector<double> r_ox(3);
    //std::vector<double> r_oy(3);
    r_ox = chain[origin]->get_position();
    //r_oy = chain[origin_link]->get_position();

    std::vector<double> u(3);
    u = vector_subtraction(r_ox, running_centre);// our rotation axis
    //normalize(u);
    //w = vector_subtraction(r_oy, running_centre);

    //double theta{ 30 * atan(1) * 4 / 180 };//30 degree rotation
    double theta{ angle };

    unpack_region(s);
    std::vector<std::vector<double>> rotated_positions(s.size());
    std::vector<double> temp(3);
    for (int i{ 0 }; i < s.size();i++) {
        temp = chain[s[i]]->get_position();
        rotated_positions[i] = rotate_by_angle_about_general_axis(temp, u, theta, com);
        //rotated_positions[i] = vector_addition(com, rotated_positions[i]);
        //chain[s[i]]->change_position(rotated_positions[i]);
        regrown_positions[s[i]] = rotated_positions[i];

    }

    // should test that with our calculated v we do indeed generate the desired helix.
    //std::vector<std::vector<double>> s_x, s_y;
    //generate(s.size() / 2, new_v, u, r_ox, s_x, s_y);


    std::vector<std::vector<int>> growth_limits;
    swivel_growth_limits(dh, growth_limits);

    grow_limits(growth_limits, alpha);

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
    std::cout << "centre of mass vector " << std::endl;
    print_1d_doub_vec(r_com);
    return r_com;
}

void polymer::u_rotate_double_helix(helix_struct* dh, double angle) {

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
    p2 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);

    regrown_positions[m1] = p1, regrown_positions[m2] = p2;
    //chain[m1]->change_position(p1), chain[m2]->change_position(p2);
    for (int i = 1; i < n; i++)
    {
        running_centre = vector_addition(running_centre, temp);

        m1 = m1 + gamma;
        m2 = m2 - gamma;

        p1 = chain[m1]->get_position();
        p2 = chain[m2]->get_position();

        p1 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);
        p2 = rotate_by_angle_about_general_axis(p1, v, angle, running_centre);

        regrown_positions[m1] = p1, regrown_positions[m2] = p2;
        //chain[m1]->change_position(p1), chain[m2]->change_position(p2);


    }

    std::vector<std::vector<int>> growth_limits;
    swivel_growth_limits(dh, growth_limits);

    grow_limits(growth_limits, alpha);

}

bool polymer::reject_o_translation(helix_struct* dh, std::vector<double>& translation) {
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

void polymer::origin_translation(helix_struct* dh, std::vector<double> t) {
    std::vector<int> s{ dh->get_monomers() };
    int alpha{ dh->get_alpha() };
    unpack_region(s);

    std::vector<double> temp(3);
    for (auto i : s) {
        temp = chain[i]->get_position();
        regrown_positions[i] = vector_addition(temp, t);
    }
    std::vector<std::vector<int>> growth_limits;
    swivel_growth_limits(dh, growth_limits);

    grow_limits(growth_limits, alpha);


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