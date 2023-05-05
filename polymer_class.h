#pragma once
#pragma once

#include "monomer_class.h"
#include "helix_struct.h"
#include "unbound_section_class.h"
#include "polymer_generation.h"
//#include "rosenbluth_sampling.h"
#include "helix_generation.h"
#include "math_functions.h"

//#include "structure_class.h"
//#include "math_functions.h"



class polymer {
protected:
    monomers chain;
    helices helix_list;
    std::vector<std::vector<int>> search_results;
    std::vector<int> extendable_structures;

    std::vector<int> zipped_structures;

    std::vector<std::vector<double>> excluded_volume;

    unbound_section* hairpin_section;
    unbound_section* free_section;

private:
    int N{ 0 };

public:

    //constructors
    polymer() = default;
    polymer(int num);

    polymer(std::vector<double> beginning_position, std::vector<double> end_position);
    polymer(std::vector<std::vector<double>> generated_positions);
    //destructor
    ~polymer();

    //general functions of the polymer and given monomer in the chain
    int chain_length();
    void get_monomer_position(int id);
    void add_monomer(std::vector<double> position, std::string base_type);
    std::vector<std::vector<double>> generate_vector_of_monomer_positions();
    void print_monomer_positions();
    void print_bases();
    void print_bases1();
    void output_for_ovito4(std::string filename="");
    void output_for_ovito5();
    void output_for_ovito6(std::string filename);
    //operator functions
    monomer* operator[](int i);

    //compatibility
    //bool base_pairing_rules(char b1, char b2);//could go to helix header
    //std::vector<int> position_of_region(int i, int r);
    //bool regions_compatible(std::string r1, std::string r2);//could go to helix header

    std::vector<std::vector<int>> structure_search(int n=3);
    std::vector<int> sample_structure(std::vector<std::vector<std::vector<int>>> potential_structures);

    bool compatible_bases(monomer* monomer_i, monomer* monomer_j);
    ////////////////////////////////////////////////////////////////////////////////////////

    void unpack_region(std::vector<int>& reaction_region);

    void pack_region(std::vector<int>& reaction_region);

    void get_linked_monomers(std::vector<int>& links);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    bool overlapping_monomers(std::vector<int> &reaction_region);
    bool reject_link(std::vector<int> &link_region, int &alpha, int&beta);
    void link();
    void unlink();

    void structure_extension(std::vector<int> &s, int &index);

    void update_large_struct_list();

    void update_extensible_structures();


    bool zip_structure_overlap(std::vector<int> s_z, int side);



    void zip_growth_limits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits);
    bool reject_zip(std::vector<int>& s, std::vector < std::vector<double>> new_positions, int side);
    void zip(bool &success);
    void sample_zip_move();
    void sample_unzip_region(std::vector<int>& region_to_unzip, int side);


    void unzip(bool& success);
    //void reject_unzip(std::vector<int>& s, int alpha, int beta);

    void unzip_growth_linits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits);

    void sample_zip_region(std::vector<int>& region_to_zip, std::vector<int>&extension, int side);

    void neighbouring_linkers();

    void grow_limits(std::vector<std::vector<int>> &limits, int alpha);

    ////////// for the simplistic hairpin simulation ////////////////////////
    bool linked();
    void link1(std::vector<int>& link_region, int alpha, int beta, int ss_index);
    void sample_link_region(std::vector<int>& link, int& alpha, int& beta, int& ss_index);
    bool reject_link1(std::vector<int>& link_region, int alpha, int beta);
    double acceptance_link0(bool link_or_unlink);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    helix_struct* sample_double_helix();
    void add_bp_to_helix(helix_struct* double_helix, std::vector<int> extension);
    std::vector<std::vector<double>> add_bp_to_helix(helix_struct* double_helix, int side);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    void update_excluded_volume(std::vector<std::vector<int>>& growth_limits, std::vector<int> helix = {});
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    void sample_swivel(helix_struct* double_helix, std::vector<double>& u_n, std::vector<double>& v_n, std::vector<double>& origin);
    void reject_swivel(helix_struct* double_helix, std::vector<double>& new_vector, int u_v_or_o);
    void swivel();
};

