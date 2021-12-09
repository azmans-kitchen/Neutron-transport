#Method of Characteristics

This folder contains scripts I wrote during my M.Sc. thesis. 

It contains the functions for solving neutron transport using Method of Characteristics(MOC)
	-> in one and in two dimensions (1D & 2D), 
	-> with one-group approximation and with multigroup approximation(1g & mg).


All the scripts were written in python language using the IDE of Spyder in Windows 10. The dependencies are numpy and multiprocessing libraries.

The 2D functions are written for a rectangular grid. To solve a problem:
	-> Describe the problem in a script
	-> Save it as "input.py".
	-> Copy the relevant "solver" in the same directory as input.
	-> Copy "raytracer2.py" in the same directory as input if the problem is in two dimensions.
	-> Run the script "solver" in that directory on Spyder. 

The eigenvalues will be printed for each source iteration. The final flux values will be stored as "phi" in the console.

To understand how to write "input" for a problem, please refer to the "input" files of the benchmarks in:

	->...\MOC\benchmark_examples\1_one_group_one_dimension\Homogenous_Slab
	->...\MOC\benchmark_examples\2_multi_group_one_dimension\URR-3-0-IN
	->...\MOC\benchmark_examples\3_one_group_two_dimension\one_group_eigenvalue_problem

To achieve optimum result with parallelization, set the number of azimuthal divisions as an integer multiple of 4*cores. 

If a parallel execution is not completed, kill the individual python processes in th task manager

Please see the reference for each of the benchmarks provided there. 

For any queries, please contact me on azmanrafee@gmail.com



