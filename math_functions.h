#pragma once
#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>
#include <cfloat>

static const double pi{ atan(1) * 4 };

double ran2(std::int64_t* idum);
double rand2(double min_value, double max_value);
//VECTOR FUNCTIONS
std::vector<double> vector_addition(std::vector<double> &vector1, std::vector<double> &vector2);
std::vector<double> vector_subtraction(std::vector<double> &vector1, std::vector<double> &vector2);
double vector_modulus(std::vector<double> &vector1);
double dot_product(std::vector<double> &vector1, std::vector<double> &vector2);
double angle_between_vectors(std::vector<double> &vector1, std::vector<double> &vector2);
std::vector<double> normalize(std::vector<double> &unnormalized_vector);
std::vector<double> multiplication_by_scalar(double constant, std::vector<double> &vector);
double dist_2_points3d(std::vector<double> &v1, std::vector<double> &v2);
std::vector<double> rotate_by_angle_about_general_axis(std::vector<double> point, std::vector<double> axis, double theta, std::vector<double> origin);

std::vector<double> cross_product(std::vector<double> &a, std::vector<double> &b);

std::vector<double> sample_unit_perp_vector(std::vector<double> u);

//OTHER
double factorial(int n);
double choose(int n, int  k);
double sum_of_elements(std::vector<double> v);
int sum_of_elements(std::vector<int> v);
std::vector<double> generate_linear_array(double min, double max, int numBins);

double product_of_elements(std::vector < double> &v);
double average_of_elements(std::vector<double>& v);
void print_1d_int_vec(std::vector<int> vec);
void print_1d_doub_vec(std::vector<double> vec);
void print_2d_int_vec(std::vector<std::vector<int>> vec);
void print_2d_doub_vec(std::vector<std::vector<double>> vec);
void double_vector_to_txt(const std::string& filename, const std::vector<double>& vector);
void int_vector_to_txt(const std::string& filename, const std::vector<int>& vector);

