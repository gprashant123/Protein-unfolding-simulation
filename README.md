# Protein-unfolding-simulation
Using Monte Carlo Simulations to unfold a 16-mer protein in a 2D lattice
## ProteinFold.m
This function initializes the native state of the protein and finds the energy of the current conformation, having arguments 
e (Energy of a single non-covalent interaction) and kT(product of k and temperature). Further, it simulates the unfolding of proteins using
Monte Carlo Method. A plot is displayed after each iteration to show the
current state of the protein. This function makes use of two functions
Energy() to calculate the energy of the native state protein and EnergyNew() to calculate the energy of the current state 

## Energy.m 

Function to calculate the energy of the native state of the protein.

Arguments:

e - Energy of a single non-covalent interaction

x and y - coordinates of the native state protein

returns:

E - The calculated total energy of the native state

new - A cell array which consists of the coordinates

amino - A 16x3 cell array which keeps track of the native state
interactions

## EnergyNew.m

Function to calculate the energy of the current conformation of the protein.

Arguments:

e - Energy of a single non-covalent interaction

x and y - coordinates of the native state protein

amino - A 16x3 cell array which keeps track of the native state
interactions 

Returns:

E - The calculated total energy of the current conformation

A detailed report along with the basic documentation of the code is given the PDF file.
