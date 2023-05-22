#pragma once
#pragma once

#include "monomer_class.h"
#include "helix_struct.h"
#include "rosenbluth_growth_class.h"
#include "polymer_generation.h"
#include "helix_generation.h"
#include "math_functions.h"




class polymer {
protected:
    monomers chain;
    helices helix_list;
    std::vector<std::vector<int>> search_results;
    std::vector<int> extendable_structures;

    std::vector<int> zipped_structures;

    std::vector<std::vector<double>> excluded_volume;

    rosenbluth_growth* hairpin_section;
    rosenbluth_growth* free_section;

    // data for acceptance probability calculation.
    std::vector<std::vector<int>> regrowth_lims;
    std::vector<std::vector < double >> regrown_positions;
    rosenbluth_growth* accepted_weights;
    rosenbluth_growth* proposed_weights;

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


    std::vector<std::vector<int>> structure_search(int n=3);
    std::vector<int> sample_structure(std::vector<std::vector<std::vector<int>>> potential_structures);

    bool compatible_bases(monomer* monomer_i, monomer* monomer_j);
    ////////////////////////////////////////////////////////////////////////////////////////

    void unpack_region(std::vector<int>& reaction_region);

    void pack_region(std::vector<int>& reaction_region);

    void get_linked_monomers(std::vector<int>& links);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    bool overlapping_monomers(std::vector<int> &reaction_region);
    void structure_extension(std::vector<int> &s, int &index);

    void update_large_struct_list();

    void update_extensible_structures();

    bool zip_structure_overlap(std::vector<int> s_z, int side);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void unlink();
    void unlink_update(int s_index); // if move is accepted
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool linked();
    void sample_link_region(std::vector<int>& link, int& alpha, int& beta, int& ss_index);
    bool reject_link1(std::vector<int>& link_region, int alpha, int beta);
    void link(std::vector<int>& link_region, int alpha, int beta, int ss_index);
    void link_update(int ss_index, int a, int b, std::vector<double>&v, std::vector<std::vector<double>> &rcentres);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void zip_growth_limits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits);
    bool reject_zip(std::vector<int>& s, std::vector < std::vector<double>> new_positions, int side);
    void zip(bool &success);
    void sample_zip_move();
    void sample_zip_region(std::vector<int>& region_to_zip, std::vector<int>& extension, int side);
    void zip_update();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void sample_unzip_region(std::vector<int>& region_to_unzip, int side);
    void unzip(bool& success);
    void unzip_growth_linits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits);
    void unzip_update();

    ////////// for the simplistic hairpin simulation ////////////////////////
    double acceptance_link0(bool link_or_unlink);
    double regeneration_probabilities();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    helix_struct* sample_double_helix();
    std::vector<std::vector<double>> add_bp_to_helix(helix_struct* double_helix, int side);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    void update_excluded_volume(std::vector<std::vector<int>>& growth_limits, std::vector<int> helix = {});
    void update_positions();
    void neighbouring_linkers();
    void grow_limits(std::vector<std::vector<int>>& limits, int alpha);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    void swivel(bool success);
    bool reject_v_swivel(helix_struct* double_helix, double angle);
    bool reject_u_swivel(helix_struct* double_helix, double angle);
    bool reject_o_translation(helix_struct* double_helix, std::vector<double> &translation);
    void v_rotate_double_helix(helix_struct* dh, double angle);
    void u_rotate_double_helix(helix_struct* dh, double angle);
    void origin_translation(helix_struct* dh, std::vector<double> t);
    void swivel_growth_limits(helix_struct* dh, std::vector<std::vector<int>>& limits);
    std::vector<double> centre_of_mass(std::vector<int> structure);
};

