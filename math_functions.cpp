#include "math_functions.h"
//#include "polymer_generation.h"


std::vector<double> vector_addition(std::vector<double> &vector1, std::vector<double> &vector2) {

    return { vector1[0] + vector2[0],vector1[1] + vector2[1],vector1[2] + vector2[2] };
}
std::vector<double> vector_subtraction(std::vector<double> &vector1, std::vector<double> &vector2) {
    return { vector1[0] - vector2[0],vector1[1] - vector2[1],vector1[2] - vector2[2] };

}

double vector_modulus(std::vector<double> &vector1) {
    return sqrt(pow(abs(vector1[0]), 2) + pow(abs(vector1[1]), 2) + pow(abs(vector1[2]), 2));
}

double dot_product(std::vector<double> &vector1, std::vector<double> &vector2) {
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}
double angle_between_vectors(std::vector<double> &vector1, std::vector<double> &vector2) {

    return acos(dot_product(vector1, vector2) / (vector_modulus(vector1) * vector_modulus(vector2)));
}


int quadrant_of_a_2d_vector(std::vector<double> two_d_vector) {
    if (two_d_vector[0] == 0.0 || two_d_vector[1] == 0) { return 0; }
    if (two_d_vector[0] > 0.0 && two_d_vector[1] > 0.0) { return 1; }
    if (two_d_vector[0] > 0.0 && two_d_vector[1] < 0.0) { return 4; }
    if (two_d_vector[0] < 0.0 && two_d_vector[1] > 0.0) { return 2; }
    if (two_d_vector[0] < 0.0 && two_d_vector[1] < 0.0) { return 3; }
}


std::vector<double> normalize(std::vector<double> &unnormalized_vector) {
    std::vector<double> unit_vector(3);
    //for (auto i : unnormalized_vector) {
    //    i = i / vector_modulus(unnormalized_vector);
    //}
    for (int i{ 0 }; i < unnormalized_vector.size(); i++) {
        unit_vector[i]=unnormalized_vector[i] / vector_modulus(unnormalized_vector);
    }
    unnormalized_vector = unit_vector;
    return unit_vector;

}
std::vector<double> multiplication_by_scalar(double constant, std::vector<double> &vector) {
    return { constant * vector[0],constant * vector[1],constant * vector[2] };
}
double dist_2_points3d(std::vector<double> &v1, std::vector<double> &v2)
{
    return sqrt(pow(v1[0]-v2[0],2)+ pow(v1[1] - v2[1], 2) + pow(v1[2] - v2[2], 2));
}
std::vector<double> cross_product(std::vector<double> &a, std::vector<double> &b) {
    return { a[1] * b[2] - a[2] * b[1],a[2] * b[0] - a[0] * b[2],a[0] * b[1] - a[1] * b[0] };
}

// probably more efficient way to do it.
std::vector<double> sample_unit_perp_vector(std::vector<double> u) {
    double mod{2.0}, cos{2.0};
    while(mod>1 or cos>0.99){
        std::vector<double> axis{ rand2(0,1),rand2(0,1),rand2(0,1) };
        mod=vector_modulus(axis);
        axis=normalize(axis);
        u=normalize(u);
        cos=dot_product(axis,u);
    }
    std::vector<double> v{ cross_product(axis,u) };
    return normalize(v);
}

#include <cmath>

// Rotate a point or vector by an angle theta about an axis u.
//chatGPT.
std::vector<double> rotate_by_angle_about_general_axis(std::vector<double> point, std::vector<double> axis, double theta, std::vector<double> origin) {

    // Translate the point to the origin.
    std::vector<double> translatedPoint(3);
    for (int i = 0; i < 3; i++) {
        translatedPoint[i] = point[i] - origin[i];
    }
    // Normalize the axis.
    double axisLength = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    std::vector<double> normalizedAxis = { axis[0] / axisLength, axis[1] / axisLength, axis[2] / axisLength };

    // Calculate the sin and cos of theta.
    double sinTheta = sin(theta);
    double cosTheta = cos(theta);

    // Calculate the components of the rotation matrix.
    double ux = normalizedAxis[0];
    double uy = normalizedAxis[1];
    double uz = normalizedAxis[2];
    double oneMinusCosTheta = 1 - cosTheta;

    // Calculate the rotation matrix.
    double rotation[3][3] = { {cosTheta + ux * ux * oneMinusCosTheta, ux * uy * oneMinusCosTheta - uz * sinTheta, ux * uz * oneMinusCosTheta + uy * sinTheta},
                             {uy * ux * oneMinusCosTheta + uz * sinTheta, cosTheta + uy * uy * oneMinusCosTheta, uy * uz * oneMinusCosTheta - ux * sinTheta},
                             {uz * ux * oneMinusCosTheta - uy * sinTheta, uz * uy * oneMinusCosTheta + ux * sinTheta, cosTheta + uz * uz * oneMinusCosTheta} };

    // Apply the rotation to the point.
    std::vector<double> rotatedPoint(3);
    for (int i = 0; i < 3; i++) {
        rotatedPoint[i] = rotation[i][0] * translatedPoint[0] + rotation[i][1] * translatedPoint[1] + rotation[i][2] * translatedPoint[2];
    }
    rotatedPoint = vector_addition(rotatedPoint, origin);
    return rotatedPoint;
}

//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//

double factorial(int n) {
    if (n == 0) return 1;
    //else if (n > 50) {
    //    std::cout << "Overload from too big a factorial";
    //    exit(0);
    //}
    else if (n < 0) {
        std::cout << "Negative integer factorial called." << std::endl;
        exit(0);
    }
    else {
        double summand{ 1 };
        for(int i = 1; i <= n; ++i)
            summand *= i;
        return summand;
    }
}
double choose(int n, int  k) {
    if (k == 0) return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

double sum_of_elements(std::vector<double> v)
{
    double sum{0};
    for (auto& n : v)
        sum += n;
    return sum;
}

double product_of_elements(std::vector<double>& v) {
    double prod{ 1 };

    for (auto& n : v) {
        prod *= n;
    }
    return prod;
}
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//

double factorial2(int n) {
    if (n == 0) return 1;
    //else if (n > 50) {
    //    std::cout << "Overload from too big a factorial";
    //    exit(0);
    //}
    else if (n < 0) {
        std::cout << "Negative integer factorial called." << std::endl;
        exit(0);
    }
    else {
        double summand{ 1 };
        for (int i = 1; i <= n; ++i)
            summand *= i * pow(10, -300);
        summand = summand;
        int p{ 300 };
        return summand;
    }
}

//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ran2.cpp

// C headers

// C++ headers

// Athena++ headers
//#include "../athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn double ran2(std::int64_t *idum)
//! \brief  Extracted from the Numerical Recipes in C (version 2) code. Modified
//!  to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
//!
//! Long period (> 2 x 10^{18}) random number generator of L'Ecuyer with Bays-Durham
//! shuffle and added safeguards.  Returns a uniform random deviate between 0.0 and 1.0
//! (exclusive of the endpoint values).  Call with idum = a negative integer to
//! initialize; thereafter, do not alter idum between successive deviates in a sequence.
//! RNMX should appriximate the largest floating-point value that is less than 1.

#define IMR1 2147483563
#define IMR2 2147483399
#define AM (1.0/IMR1)
#define IMM1 (IMR1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

double ran2(std::int64_t* idum) {
    int j;
    std::int64_t k;
    static std::int64_t idum2 = 123456789;
    static std::int64_t iy = 0;
    static std::int64_t iv[NTAB];
#pragma omp threadprivate(iy,iv,idum2)
    double temp;

    if (*idum <= 0) { // Initialize
        if (-(*idum) < 1)
            *idum = 1; // Be sure to prevent idum = 0
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) { // Load the shuffle table (after 8 warm-ups)
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0) *idum += IMR1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;                 // Start here when not initializing
    *idum = IA1 * (*idum - k * IQ1) - k * IR1; // Compute idum=(IA1*idum) % IMR1 without
    if (*idum < 0) *idum += IMR1;   // overflows by Schrage's method
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2; // Compute idum2=(IA2*idum) % IMR2 likewise
    if (idum2 < 0) idum2 += IMR2;
    j = static_cast<int>(iy / NDIV);              // Will be in the range 0...NTAB-1
    iy = iv[j] - idum2;                // Here idum is shuffled, idum and idum2
    iv[j] = *idum;                 // are combined to generate output
    if (iy < 1)
        iy += IMM1;

    if ((temp = AM * iy) > RNMX)
        return RNMX; // No endpoint values
    else
        return temp;
}

#undef IMR1
#undef IMR2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


double rand2(double min_value, double max_value) {

    std::int64_t seed{ 5 };
    return double(((max_value - min_value) * ran2(&seed))) + min_value;
}


///////////**********************************************************************//////////////////////
///////////**********************************************************************//////////////////////
///////////**********************************************************************//////////////////////
void print_1d_int_vec(std::vector<int> vec)
{
    for (auto i : vec) {
        std::cout << i << " ";
    }std::cout << std::endl;

}

void print_1d_doub_vec(std::vector<double> vec)
{
    for (auto i : vec) {
        std::cout << i << " ";
    }std::cout << std::endl;
}

void print_2d_int_vec(std::vector<std::vector<int>> vec)
{
    for (auto i : vec) {
        for (auto j : i) {
            std::cout << j << " ";
        }std::cout << std::endl;
    }std::cout << std::endl;
}

void print_2d_doub_vec(std::vector<std::vector<double>> vec)
{
    for (auto i : vec) {
        for (auto j : i) {
            std::cout << j << " ";
        }std::cout << std::endl;
    }std::cout << std::endl;

}
///////////**********************************************************************//////////////////////
///////////**********************************************************************//////////////////////
///////////**********************************************************************//////////////////////
