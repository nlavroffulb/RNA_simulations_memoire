#include "helix_generation.h"
#include "math_functions.h"
#include "polymer_generation.h"
//#include "rosenbluth_growth_class.h"

std::string random_base() {
    int rn = rand() % 4;
    std::string base[4] = { "C","U","A","G" };

    return base[rn];
}



bool base_pairing_rules(char b1, char b2)
{
    if (b1 == 'A' && b2 == 'U') {
        return true;
    }
    else if (b1 == 'U' && b2 == 'A') {
        return true;
    }
    else if (b1 == 'C' && b2 == 'G') {
        return true;
    }
    else if (b1 == 'G' && b2 == 'C') {
        return true;
    }
    else return false;
}

bool regions_compatible(std::string r1, std::string r2)
{
    int counter{ 0 };
    //std::cout << "length of string " << r1.length() << std::endl;
    for (int i{ 0 }; i < r1.length(); i++) {
        base_pairing_rules(r1.at(i), r2.at(r1.length() - i - 1)) ? counter++ : counter = 0;
    }
    if (counter == r1.size()) {
        //std::cout << "Compatible" << std::endl;
        return true;
    }
    else {
        return false;

    }
}

std::vector<int> position_of_region(int i, int r)
{
    std::vector<int> reg;
    for (int j{ 0 }; j < r; j++) {
        reg.push_back(i + j);
    }
    return reg;
}








//major groove
//const double Mg = 0.47;
const double Mg = 50 * (0.47 / 1.08);

//minor groove
//const double mg = 1.08;
const double mg = 50;

//diameter
//const double dm = 2.5;
const double dm = 50;
//const double theta = 0.75 * atan(1) * 4;
//const double phi = 33 * (atan(1) * 4) / 180;
//const double radius{ 0.5 };

const double theta = 0.75 * atan(1) * 4;
const double phi = 33 * (atan(1) * 4) / 180;
const double radius{ 0.25 };

const double pitch{ 0.28 };

double the = acos((2 * dm * dm / 4 - mg * mg) / (2 * dm * dm / 4.));
//print('the=', the, the * 180 / np.pi)
void sample_helix_vectors(std::vector<double>& u, std::vector<double>& v) {
    sample_jump_direction(v, 1);

    u = sample_unit_perp_vector(v);

}

//returns the distance between bonded monomers in the helix ie the diameter. this is useful for calculating
//overstretch condition later. 
double helix_separation() {
    std::vector<double> u, v, x{ rand2(0,1),rand2(0,1),rand2(0,1) };
    sample_helix_vectors(u, v);
    std::vector<double> temp{ multiplication_by_scalar(radius,u) };
    std::vector<double> running_centre{ vector_subtraction(x,temp) };

    std::vector<double> w{ rotate_helix(u,v,theta) };

    temp =  multiplication_by_scalar(radius, w);
    std::vector<double> y = vector_addition(running_centre, temp);

    return dist_2_points3d(x, y);

}

//rotate u around v by angle theta
std::vector<double> rotate_helix(std::vector<double>& u, std::vector<double>& v, double the) {// v is the axis of rotation
    normalize(v);
    double cthe = cos(the);
    double sthe = sin(the);
    double v_dot_u = dot_product(u, v);
    //v_vec_u = np.zeros(3)
    std::vector<double> v_vec_u(3);

    //cross product
    v_vec_u[0] = v[1] * u[2] - v[2] * u[1];
    v_vec_u[1]=v[2] * u[0] - v[0] * u[2];
    v_vec_u[2]=v[0] * u[1] - v[1] * u[0];

    std::vector<double> wx;
    wx = multiplication_by_scalar(cthe, u);
    std::vector<double> temp1{ multiplication_by_scalar(sthe, v_vec_u) }, temp2{ multiplication_by_scalar((1.0 - cthe) * v_dot_u, v) };
    wx = vector_addition(wx, temp1);
    wx = vector_addition(wx, temp2);


    return wx;

}
//
//generate an A - RNA helix of n bp
//u: helix direction
//x : euclidean position of the starting base pair backbone
//v : versor joining the center of the helix to x; uand v are orthogonal




// it is imperative that u and v be normalized and that they be orthogonal. very important otherwise we can't predict overstretch condition.
void generate(int region_length, std::vector<double>& v,
    std::vector<double>& u, std::vector<double>& x,
    std::vector<std::vector<double>>& s_x, std::vector<std::vector<double>>& s_y) {

    int n{ region_length};

    //if (compatible_monomers.size() == 4) {
    //    n= compatible_monomers
    //}

    std::vector<std::vector<double>> strand_x(n);
    std::vector<std::vector<double>> strand_y(n);
    //double theta = 0.75 * atan(1)*4;
    //double phi = 33 * (atan(1) * 4) / 180;
    //double dv{ 0.28 }, radius{ 2.5 / 2 };
    //for (auto i : x) {
    //    std::cout << i << " ";
    //}
    //std::cout << std::endl;
    std::vector<double> temp{ multiplication_by_scalar(radius,u) };
    std::vector<double> running_centre{ vector_subtraction(x,temp) };

    std::vector<double> w{ rotate_helix(u,v,theta) };

    temp = multiplication_by_scalar(radius, w);
    std::vector<double> y = vector_addition(running_centre, temp);

    strand_x[0] = x;
    strand_y[n - 1] = y;
    //for (auto i : strand_y[n - 1]) {
    //    std::cout << i << " ";
    //}
    //print_1d_doub_vec(running_centre);
    //print_1d_doub_vec(w);

    int running_n{ 1 };
    std::vector<double> running_x(3);
    std::vector<double> running_y(3);
    std::vector<double> temp_u{u};
    std::vector<double> up(3), wp(3);
    while (running_n < n) {
        std::vector<double> temp{ multiplication_by_scalar(pitch, v) };

        running_centre = vector_addition(running_centre, temp);
        //running_centre = vector_addition(running_centre, multiplication_by_scalar(0.28, v));
        up = rotate_helix(temp_u, v, phi);

        wp = rotate_helix(w, v, phi);

        temp_u = up;

        w = wp;
        temp = multiplication_by_scalar(radius, temp_u);

        running_x = vector_addition(running_centre, temp);
        temp = multiplication_by_scalar(radius, w);
        running_y = vector_addition(running_centre, temp);


        //the side of the compatible region that contains the origin of the helix is filled up
        // following the order of the monomers in the chain. = s_x
        // s_y is also filled up according to the order of the chain. 
        //if a region is 4 5 6 10 11 12 then s_x = [4 5 6], s_y = [10 11 12]
        strand_x[running_n] = running_x;
        strand_y[n - running_n - 1] = running_y;
        running_n++;

    }
    s_x = strand_x;
    s_y = strand_y;



}

void generate(int region_length,
    std::vector<double>& v, std::vector<double>& u, std::vector<double>& x,
    std::vector<std::vector<double>>& s_x, std::vector<std::vector<double>>& s_y, std::vector<std::vector<double>> &running_cs) {

    int n{ region_length };
    std::vector<std::vector<double>> strand_x(n);
    std::vector<std::vector<double>> strand_y(n);
    std::vector<double> temp{ multiplication_by_scalar(radius,u) };
    std::vector<double> running_centre{ vector_subtraction(x,temp) };


    //helix->set_running_x_strand_centre(running_centre);
    std::vector < std::vector<double>> temp_rcs(2);
    temp_rcs[0] = running_centre;

    std::vector<double> w{ rotate_helix(u,v,theta) };

    temp = multiplication_by_scalar(radius, w);
    std::vector<double> y = vector_addition(running_centre, temp);

    strand_x[0] = x;
    strand_y[n - 1] = y;
    //for (auto i : strand_y[n - 1]) {
    //    std::cout << i << " ";
    //}
    //print_1d_doub_vec(running_centre);
    //print_1d_doub_vec(w);

    int running_n{ 1 };
    std::vector<double> running_x(3);
    std::vector<double> running_y(3);
    std::vector<double> temp_u{ u };
    std::vector<double> up(3), wp(3);
    while (running_n < n) {
        std::vector<double> temp{ multiplication_by_scalar(pitch, v) };

        running_centre = vector_addition(running_centre, temp);
        //running_centre = vector_addition(running_centre, multiplication_by_scalar(0.28, v));
        up = rotate_helix(temp_u, v, phi);

        wp = rotate_helix(w, v, phi);

        temp_u = up;

        w = wp;
        temp = multiplication_by_scalar(radius, temp_u);

        running_x = vector_addition(running_centre, temp);
        temp = multiplication_by_scalar(radius, w);
        running_y = vector_addition(running_centre, temp);


        //the side of the compatible region that contains the origin of the helix is filled up
        // following the order of the monomers in the chain. = s_x
        // s_y is also filled up according to the order of the chain. 
        //if a region is 4 5 6 10 11 12 then s_x = [4 5 6], s_y = [10 11 12]
        strand_x[running_n] = running_x;
        strand_y[n - running_n - 1] = running_y;
        running_n++;

    }
    //helix->set_running_y_strand_centre(running_centre);
    temp_rcs[1] = running_centre;
    running_cs = temp_rcs;


    s_x = strand_x;
    s_y = strand_y;



}





double side_length(int n) {

    std::vector<double> u{ 0,0,1 }, v{ 1,0,0 }, o{ 0,0,0 };
    std::vector<double> temp{ multiplication_by_scalar(radius,u) };
    std::vector<double> running_centre{ vector_subtraction(o,temp) };
    temp = multiplication_by_scalar(pitch*(n - 1), v);
    running_centre = vector_addition(running_centre, temp);
    std::vector<double> up = rotate_helix(u, v, (n-1)*phi);
    temp = multiplication_by_scalar(radius, up);
    std::vector<double> running_x{ vector_addition(running_centre,temp) };

    return dist_2_points3d(o, running_x);

}

double sideways_length(int n)
{
    std::vector<std::vector<double>> s_x(n), s_y(n);
    std::vector<double> u{0,0,1}, v{0,1,0}, o{0,0,0};
    generate(n, v, u, o, s_x, s_y);

    return dist_2_points3d(s_x[0], s_y[0]);

    //return sqrt(pow(side_length(n), 2) + pow(helix_separation(), 2));
}

