# pipette_opt
Latest change: 1 July 2020 - IDDFS now works (though slowly), common auxiliary funcrions between methods taken to a separate file auxil.py

The tsp_method program implements the TSP-based algorithm outlined in the 'PossibleSolution' document. The method used to solve the TSP is 2-opt of the Nearest Neighbour algorithm

The statespace_methods program implements various algorithms working with the state-space representation of the problem (see the 'ProblemRepresentation' document). Currently, it includes an iddfs solver of invariable depth 1

The input_generator program contains fucntions that create random inputs used for algorithm testing
