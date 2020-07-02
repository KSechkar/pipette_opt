# pipette_opt
Latest change: 2 July 2020 - an alternative cost function introduced: thanks to it, an adjusted IDDFS works faster, so a depth 2 IDDFS is now possible.
The heurstic reorderings are now transferred to tsp_reorder.py from tsp_method.py of the branch 'heuristic_reorderings' (the latter is now closed).

The tsp_method program implements the TSP-based algorithm outlined in the 'PossibleSolution' document.
The method used to solve the TSP is 2-opt of the Nearest Neighbour algorithm.

The statespace_methods program implements various algorithms working with the state-space representation of the problem (see the 'ProblemRepresentation' document).
Currently, it contains two kinds of iddfs solver.

The input_generator program contains fucntions that create random inputs used for algorithm testing.

The auxil program contains auxiliary functions used by BOTH methods
Currently, it has a display function, as well as route/single operation cost functions of two different kinds.
