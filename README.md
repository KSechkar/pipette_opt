# pipette_opt

Latest change: 20 August 2020

The algorithms do work in a situation where every vector has its own required volume, while the pipette capacity is limited.
Any munber of vector types is supported

The hub_spoke program implemets the hub-and-spoke method for solving the problem.

The tsp_method program implements the LP-based algorithm, whose non-capacitated version outlined in the 'PossibleSolution' document.
The linear programming problem that is currently solved instead of the TSP is given in detail in the 'CapacityAndTSPmethod' document.

The tsp_reorder program contains the reordering algorithms that can be performed to improve the preformance of tsp_method.
It has simple and state-space reorderings, explained in the 'Reorderings' document.

The tsp_lp_solver program implements linear programming solvers needed for the tsp-based methods.
tsp_lp_gurobi solves the TSP (for non-capacitated problem).
lp_cap uses gurobi to solve the linear programming problem needed for the capacity-conscious version of the 'TSP method'.

The statespace_methods program implements various algorithms working with the state-space representation of the problem (see the 'ProblemRepresentation' document).
Currently, it contains two kinds of iddfs solver and greedy solver. Monte-Carlo and A* solver are not supported any more.

The input_generator program contains fucntions that create random inputs used for algorithm testing.

The auxil program contains auxiliary functions used by BOTH methods.
Currently, it has a display function, pipette capacity calculator, as well as route and single operation cost functions.
There are also pipette capacity calculator functions. To ignore capacity, give caps=None.

The pipette_opt program contains API-like functions that allow to integrate the algorithms into the dna_assembler and DNABOT assembly pipelines.
(WORK IN PROGRESS)

***
Development purposes only

howto_manual.txt outlines how to launch and test solvers

The test_inp program creates test input .csv files

The testing program runs on .csv inputs and writes down results
