# C-code-Molecular-Dynamics-Simulation
 Molecular dynamics simulation using C code
 
 Molecular dynamic using 3D Lennard - Jones potential. Using periodic conditions, gaussian distribution for temperature and Verlet algorithm we perform a minimization of potential energy. Then, we include Langevin thermostat to get forces between particles.

 To run this code, you must clone this repository and compile these into an executable:

 gcc -Wall -O3 -o molecular_dynamic.e molecular_dynamic.c -lm

 Once you have the executable, i can be run from command line:

 ./molecular_dynamic.e

 NOTE: if oyu change a parameter in molecular_dynamic.c you must recompile the code.
