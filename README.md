# RNA_simulations_memoire

The main class is the polymer class. 
The polymer class has double helix structures which are described by an object called helix_struct. It has a header and 
cpp for function definitions. 
The polymer class also has unbound sections which are described by an object called unbound_section. This class is still in progress. Ideally
it would store the Rosenbluth weights and energies to be used in the calculation of the acceptance probability.
Thirdly, polymer contains monomer objects which are a quite basic class, just for storing the position and a couple of other properties
associated to each monomer.

Then we have polymer generation, helix generation header and cpp files, contents follow their name.
Math functions header and cpp, contents also follow the name. Mostly vector functions.
