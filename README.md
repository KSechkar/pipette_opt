# pipette_opt
Latest change: 30 June 2020 - json_reader branch merged; introduced a primitive iddfs method

The tsp_method program implements the TSP-based algorithm outlined in the 'PossibleSolution' document. The method used to solve the TSP is 2-opt of the Nearest Neighbour algorithm

The statespace_methods program implements various algorithms working with the state-space representation of the problem (see the 'ProblemRepresentation' document). Currently, it includes an iddfs solver of invariable depth 1
---
json_reader branch

Instead of a sample input list w, read a .json file and transcribe it into the list susbsets of type Ss
A dictionary reagdic maps reagent names such as J23100 into e.g. p1
---
