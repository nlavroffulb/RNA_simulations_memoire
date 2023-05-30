// the last dance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include "math_functions.h"
#include "helix_generation.h"
#include "simulations.h"
//#include "acceptance.h"
#include<time.h>
#include <fstream>
#include <chrono>
#include <algorithm>

//#include "ran2.h"

// average time for 2 consecutive link moves: 1937 ms (almost 2 seconds yikes). 4/4 15:00. need to get this down a lot.
// 1170ms 5/4 20:00 KEEP GOING !!!
// 830ms 5/4 20:30 RAHHHHHHHHHH!!!

// average time for 3 consecutive link moves 27/4: 700ms. quite good.

// time for 2 links and 1 unlink
// 1356ms 5/4 21:00
// 
// 
// 
// 
//TO DO.  
////////////////////////////////////////////////////////////////////////////////////////

// excluded volume of helixes. could still be done.
// 
// want to make an off switch for Rosenbluth sampling. I would guess the simplest thing to do would be to set k = 1.
// 
// how does swivel work with zipped structures? before we assumed that the centre of mass was in the plane of 
// one of the base pairs however if we have an even number of base pairs e.g 4 this will not be the case, so is 
// it working correctly? probably not.
// 
// getting a subscript vector out of range occasionally. still quite a few bugs to correct and generally 
// the acceptance rules are not clean.
// 
// can make metropolis structure more efficient. e,g only calculate search results, extensible structures when 
// a certain move has been accepted. after a series of rejections it should not be re calcualted.
// 
// might need to include normalization factor for yamakawa probability in the acceptance rules.
// 
// change zip so that we always sample the side of zip from which we can actually extend the helix. need to modify
// the acceptance rule potentially.
// 
// global constants header file.
// 
// acceptance rate is often nan. need to fix. 
// 
// do we always reject links if there are no structures left to try?
// 
// are we protected against structures which contain one of the terminal endpoints? need to check for all moves. 
// 
// unlink growth limits. sometimes we have growth limits that connect e.g [38,16] & [41,38]. should just be one growth limit.
// still need to include the extra probabilities which arise for each acceptance rule for example
// the probability of picking a zipped structure. 
// 
// structure id needs to be cleaned up i think. 
// // 
// when unlinking, is there any point in sampling alpha and beta?? the number of growth limits will be the same in either case.
// 
// limiting case for zip where we zip the two tail ends of the chain. am i protected against that?
// 
//  keep doing this thing for alpha  =  1 case where we reverse the order of the vector... is this necessary? when we do constraints
// we could just consider that w,x,y,z go from 3' to 5'. don't need to reverse it. remove the alpha dependency from w,x,y,z. then 
// there would be no need to reverse them.
// 
// for metropolis hastings: 
// ------------all our MC move functions shouldn't directly change the positions of the monomers because even if a move is successful
// it could be rejected because of the acceptance probability.
// 
////////////////////////////////////////////////////////////////////////////////////////

// 



// visualization of the structures that i make: different structures should have different colours.
// the monomers they enclose as well so we can easily follow the sequence.
// -- separate files for each structure. 
// -- different colours for connecting parts




int main()
{
    srand(time(NULL));
    std::int64_t init{ rand()};
    init = -1*init;
    ran2(&init);


    
    polymer p0(70);
    double u1, u2, u3;
    int NMC{ 1000 },i{0};

    link_unlink_hairpin_sim(NMC, 50, false);





}
