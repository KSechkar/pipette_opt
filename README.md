# pipette_opt
Latest change: 20 July 2020 - state-space method now supports limited pipette capacity

The tsp_method program implements the TSP-based algorithm outlined in the 'PossibleSolution' document.
The method used to solve the TSP is 2-opt of the Nearest Neighbour algorithm.

The tsp_reorder program contains the reordering algorithms that can be performed to improve the preformance of tsp_method.
It has simple and state-space reorderings.

The tsp_lp_solver program implements a linear programming approach to solving a TSP instead of the tspy package.

The statespace_methods program implements various algorithms working with the state-space representation of the problem (see the 'ProblemRepresentation' document).
Currently, it contains two kinds of iddfs solver, greedy solver and an a* solver (does not work).

The input_generator program contains fucntions that create random inputs used for algorithm testing.

The auxil program contains auxiliary functions used by BOTH methods.
Currently, it has a display function, pipette capacity calculator, as well as route/single operation cost functions of two different kinds.
Only cost functions 'with_w' are used as the other kind is too slow - and only 'with_w' supports capacity limitations. To ignore capacity, give cap=None.

***
Development purposes only

howto_manual.txt outlines how to launch and test solvers

The test_inp program creates test input .csv files

The testing program runs on .csv inputs and writes down results
