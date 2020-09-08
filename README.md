# pipette_opt

Latest change: 8 August 2020

---
FUNCTION OF EACH PROGRAM

-- Main programs (contain the methods)
- hub_spoke.py: the hub-and-spoke method for solving the problem.

- tsp_method.py: the LP-based algorithm (the name 'tsp' instead of 'lp' is a relic of the older versions).

- statespace_methods.py: optimisation algorithms using the state-space representation of the problem (see the 'ProblemRepresentation' document).


-- Auxiliary programs (the methods need them to work)
- tsp_reorder.py: the heuristic reordering algorithms that can be performed to improve the preformance of the LP method 'Sametogether' can also be used for state-space methods.

- tsp_lp_solver: linear programming solvers needed for the lp-based methods.
	Function tsp_lp_gurobi() solves the TSP (problem with no pipette capacity).
	lp_cap uses gurobi to solve the linear programming problem needed for the pipette capacity-conscious version of the 'TSP method'.

- auxil.py: auxiliary functions used by BOTH methods - output display, pipette capacity calculator; route and single operation cost functions.


-- Adapting the algorithms
- pipette_opt.py: adapter functions for the DNA-BOT (BASIC assembly) or dna_assembler (Start-stop assembly) automated DNA assembly packages.
	! The Hub-spoke methods are NOT supported
  	Work in progress - MoClo automation package adapters.

- vis_cont.py: visulaises the results of the adapters' works by reading a .p file created by them


-- Programs that were used in testing
- input_generator.py: contains fucntions that create random inputs used for algorithm testing in the folder named 'input'.

- testing.py: was used to test the functions with random-generated Start-Stop assembly cases (sections 4.1-4.2 of the Technical Report).

---
HOW TO TEST THE PROGRAMS?
(check howto_manual.txt for more details)

-- The main() functions in all programs except for vis_cont.py are exclusively for testing.
- main() functions of hub_spoke.py, statespace_methods.py and tsp_method.py allow to call optimisation algorithms on the input defined in them.
  Follow the guidelines in these files' comments to alter the test input, if needed.

-- Contrary to what the name suggests, testing.py is generally NOT for human testing.
   It does optimisation on numerous inputs and records the average resuts. The runtime is at least 24 hours if all numbers of wells from 2 to 96 are selected, depending on the algorithms chosen for testing.

-- To use the package in conjuction with DNA-BOT (BASIC assembly) or dna_assembler (Start-stop assembly) packages:
	1. Download the assembly automation package; pipette_opt must be in the adjacent folder
	2. Change the files in the original package as insturcted by the comments at the start of pipette_opt.py
	3. Run the main program of the original package
	4. Irrespectively of whether the program is a script itself or a script generator, one or multiple .p files will appear
	5. Open the .p file using the vis_cont.py visual app. A visualisation of the well plate, allowing you to view the order in which the pipette visited the well, will apear