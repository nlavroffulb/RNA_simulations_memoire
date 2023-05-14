#pragma once
#include<vector>
#include<iostream>
const int standard_struct_size{ 3 };
class helix_struct {

private:
	std::vector<double> u;
	std::vector<double> v;
	std::vector<int> monomer_indices;
	int origin;
	bool extendible;
	//bool large{false};
	std::vector<int> extension;

	std::vector < std::vector<double>> running_centres;
	int alpha, beta;
	int gamma; // gamma is for growth direction. it is the direction for which v is defined. for alpha beta = 00 , 11, v defines growth from
	// 3' to 5' while for the other options it defines 5' to 3'. we only need to set alpha and beta for a double helix once, then we'll 
	// store the running centres of each of the ends and calculate u' and w' respectively which are the vectors joining the running centre
	// to the positions of the 2 monomers in the at an extreme.

public: 
	//constructor
	helix_struct() = default;
	helix_struct(int o,std::vector<double> &u, std::vector<double> &v, std::vector<int> &monomer_indices, int growth_direction);
	helix_struct(std::vector<int> monomers, std::vector<double>& v, std::vector < std::vector<double>> r_centres, int alpha, int beta);
	helix_struct(std::vector<int>& monomers, std::vector<double>& v, int alpha, int beta);
	//helix_struct(std::vector<double> u, std::vector<double> v, std::vector<int> monomer_indices,std::vector<int> extension);

	~helix_struct(){
		std::cout << "Helix struct destructor called." << std::endl;
	};

	void set_growth_direction();
	void set_running_x_strand_centre(std::vector<double> rxc);
	void set_running_y_strand_centre(std::vector<double> ryc);
	void set_alpha(int alpha);
	void set_beta(int beta);

	// get functions
	std::vector<double> get_u(), get_v();
	int get_origin();
	std::vector<int> get_monomers(), get_extension();
	int get_growth_direction();


	void set_extension(std::vector<int> big_struct), set_extendable(bool ex);
	void set_origin(int o);
	void change_u(std::vector<double> u_vec);
	void change_v(std::vector<double>v_vec);
	void change_u_v(std::vector<double>& u, std::vector<double>& v);
	void extend(int o, std::vector<double> &u, std::vector<double> &v);
	void extend(std::vector<int> &x10sion);

	void shorten(int o, std::vector<double> &u, std::vector<double> &v);
	void shorten(std::vector<int>& x10sion);


	bool extended();
	bool fully_extended();
	bool extendable();

	int get_beta();
	int get_alpha();
	int get_last_y();//when the structure was created and generated, 
	//it is grown with 2 strands, s_x and s_y. this function gets the last monomer grown in the y strand.
	int get_mon_before_origin();
	std::vector<double> get_rc(int index);
	void set_rc(int index, std::vector<double> new_rc);
};

typedef std::vector<helix_struct*> helices;


