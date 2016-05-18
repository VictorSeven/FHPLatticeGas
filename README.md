# FHPLatticeGas

This is a FHP Lattice Gas Cellular Automata implementation in C++. It has been done following the steps of the book by D.A. Wolf-Gladrow,  "Lattice Gas Cellular Automata and Lattice Boltzmann Models" (Springer Berlin, 2000). 

It allows you to simulate the FHP-I, FHP-II or FHP-III model, simply setting up the initial conditions for the particles. Since it has been made using multispin coding, the implementation is really efficient (it can proccess a large amount the particles in a few seconds) and useful to do fluid simulations. The code measure the mean occupation numbers, velocity and density in each node.

Also, the code is well-structured and fully commented, so it's easy to modify it for your own purposes :D

FHP-III implementation differs from the one given in the book by Wolf-Gladrow (but it gives correct results). 
