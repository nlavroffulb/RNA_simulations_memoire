#pragma once
#include<vector>
#include<iostream>
const int standard_struct_size{ 3 };// we always start with double helices of 3 base pairs.
class helix_struct {

private:
	std::vector<double> u;
	std::vector<double> v;
	std::vector<int> monomer_indices;
	int origin;
	std::vector<int> extension;

	std::vector < std::vector<double>> running_centres;// we store the running centres particularly for zip moves, it's faster
	// to just add 2 monomers to the 'top' or 'bottom' of a helix using the running centres of that helix rather than re generate
	// the whole helix.
	int alpha, beta;// we need to store alpha and beta to know in what direction the v helix vector is directed.

public: 
	//constructor
	helix_struct() = default;
	helix_struct(std::vector<int> monomers, std::vector<double>& v, std::vector < std::vector<double>> r_centres, int alpha, int beta);

	~helix_struct(){
		std::cout << "Helix struct destructor called." << std::endl;
	};

	// get functions
	std::vector<double> get_u(), get_v();
	int get_origin();
	std::vector<int> get_monomers(), get_extension();
	std::vector<double> get_rc(int index);// get running centre.
	

	void set_extension(std::vector<int> big_struct);
	void set_origin(int o);
	void change_v(std::vector<double>v_vec);
	void extend(std::vector<int> &x10sion);

	void shorten(std::vector<int>& x10sion);


	bool extended();
	bool fully_extended();
	bool extendable();

	int get_beta();
	int get_alpha();
	void set_rc(int index, std::vector<double> new_rc);
};

typedef std::vector<helix_struct*> helices;// this allows to make a vector of helix struct objects that we use in the polymer class.


