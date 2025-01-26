This folder contains files for our final project that uses up to 35 functions to estimate total energy configurations that are likely to occur for a given dataset representing a set number of particles and computational restrictions.


The files contained in this repository are:
1. 2.py file consisting of our MC simulation using the Python Standard Library (PSL) and parameters specified to run the simulation.
2. 2 .py file containing our NumPy version of MC simulation and parameters specificed to run the simulation.
3. 1 .cpp file containing the C++ version of our MC simulation.
4. 1 README.md file which you are reading.


In order to run the .py file using only the PSL, the following packages need to be imported:`Math`, `Random`


In order to run the .py file using NumPy, the following packages need to be imported: `Math`, `NumPy`


In order to run the .cpp file the following libraries need to be imported: `iostream`, `fstream`, `array`, `vector`, `utility`, `random`, `chrono`, `stdio.h`.

To perform time analyses, the terminal's `time` command line utility was used for both Python and C++. The terminal was also used to compile and run our programs.


The main goal of this project was to optimize the performance of our programs to estimate the total energy that a given system is likely to be in using the lennard jones potential energy equation (in the form of reduced units), and different computational tools that was first programmed using Python's standard library, and then refactored using Python's NumPy library and C++ 

The computational and statistical methods used in this molecular dynamics simulation was: the Monte Carlo (MC) method, Metropolis-Hastings algorithm, the boltzmann probability distribution, the idea of a simulation box to model infinite space, periodic boundaries, and the implementation of cut-off values to reduce the computational cost needed to make estimations at infinite values. 


We ran 3 simulations with three step intervals: 2500 steps, 25000 steps, and 50000 steps. Each data point was run with an average of 3 trials, totalling to 27 simulation runs between the PSL, NumPy, and C++ variations of MC simulation.


A graph of the three variations of the MC simulations was created on google sheets to visualize the performance analysis of each. 

Our results indicate that Python standard library is much slower for our program when comapred to using NumPy libraries and C++. Moving forward, it may be best to use NumPy and C++ for molecular dynamic simulations instead of the python standard library. 