# pipette_opt

Latest change: 22 July 2020 - tsp method now supports limited pipette capacity; cost function fixed (was giving incorrect results)

The branch CAPACITY introduced pipette capacity limitations to the problem.
Currently, the algorithms do work in a situation where every vector has its own required volume (and thus own capacity).
Due to a merger with branch OTHER_NUMS the algorithms also work for any number of part types.

The hub_spoke program implemets the hub-and-spoke method for solving the problem.

The tsp_method program implements the TSP-based algorithm outlined in the 'PossibleSolution' document.
The method used to solve the TSP is 2-opt of the Nearest Neighbour algorithm.

The tsp_reorder program contains the reordering algorithms that can be performed to improve the preformance of tsp_method.
It has simple and state-space reorderings.

The tsp_lp_solver program implements linear programming solvers needed for the tsp-based methods.
cvxpy-based methods are too slow to work on 96 wells (are commented right now to reduce distraction), but tsp_lp_gurobi is fast and working well.
lp_cap uses gurobi to solve the linear programming problem needed for the capacity-conscious version of the 'TSP method'.

The statespace_methods program implements various algorithms working with the state-space representation of the problem (see the 'ProblemRepresentation' document).
Currently, it contains two kinds of iddfs solver, greedy solver and an a* solver (does not work).

The input_generator program contains fucntions that create random inputs used for algorithm testing.

The auxil program contains auxiliary functions used by BOTH methods.
Currently, it has a display function, pipette capacity calculator, as well as route/single operation cost functions of two different kinds.
There are also pipette capacity calculator functions. To ignore capacity, give caps=None.
Only cost functions 'with_w' are used as the other kind is too slow - and only 'with_w' supports capacity limitations.

***
Development purposes only

howto_manual.txt outlines how to launch and test solvers

The test_inp program creates test input .csv files

The testing program runs on .csv inputs and writes down results
