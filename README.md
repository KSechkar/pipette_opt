# pipette_opt

Latest change: 8 April 2021

This repository contains the Python 3.8 code implementing the algorithms for opitmising the pipette tip uptake during the distribution of DNA parts across the well array in which DNA assembly constructs are prepared. The detailed explanation of the algorithms can be found at <link to the article and/or preprint>.

## Algotihm implementations
There are three main files, each of them implementing one of the three approaches to solving the tip consumption optimisation problem:
* statespace_methods.py - searching a tree graph of states of the system
* lp_method.py - dividing the problem into a series of Linear Programming problems, and using the [GUROBI](https://www.gurobi.com/) optimiser to solve them
* dp_method.py - dynamic programming
* hub_spoke.py - 'hub and spoke' method (OLD VERSION; CURRENTLY NOT SUPPORTED)

All of the algorithms recieve an input in an abstract format independent from the assembly standard: each DNA part is represented as a tuple (a,b), i.e. species number b in the list of parts found on position a in the assembled contructs. The dictionary _caps_  outlines how many doses of each part's solution the pipette can hold; each nested array in the  2D-list _w_ outlines the composition of a single construct to be prepared. The output is an array _fin_, where each entry stands for the addition of a given DNA part to a single construct well; it also specifies if the tip must be changed to perform this operation.

The main() functions in these files allow to run the algorithms with certain inputs to test them. By changing and (un)commenting code lines in main, the algorithm can be run on a defined test input or on a randomly-generated input of up to 96 constructs.  Depending on what line of the code is uncommented, different algorithms from the same file can be tested out. By changing the reord argument, the reordering of the operation list (preprocessing of input to improve algorithm performance) can be selected.

The files:
* auxil.py
* lp_reorder.py
* lp_solver.py

## Adaptation for different DNA assemblies
pipette_opt.py allows to enable pipette tip saving in existing automated DNA assembly protocols for Opentrons OT-2 for [Start-Stop](https://github.com/zoltuz/dna_assembler), [BASIC](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) and [MoClo](https://github.com/DAMPLAB/OT2-MoClo-Transformation-Ecoli) standards. It contains code and insturctions how to change the original protocol files. Please note that while the conversion of input into the internal abstract format and pipette tip consumption optimisation have been validated, no DNA assemblyhas been executed and tested in accordance with the obtained protocol.

For the Start-Stop and BASIC assembly, the resultant sequence of actions is saved in a .p file. By running the vis_cont.py programme and reading this .p file, the possible contaminations between consturct wells can be visualised.

## Algorithm testing
Randomly generated inputs in the internal abstract format can be saved in .csv files:
* input_generator.py contains the functions randomly generating Start-Stop assembly inputs with a given number of wells and parts of each type and saving them in .csv files
* test_inp.py for n from 2 to 96 saves a .csv file that contains 100 inputs, n constructs each

The program testing.py reads the files created as test_inp.py and finds mean numbers of tips required and algorithm runtimes, along with standard deviations, saving them as .csv entries. The program should be run from command prompt with the following arguments
* -w which algorithms to test:
	* all - all algorithms
	* statespace - all state-space algorithms
		* sssametogether - state-space with 'sametogether' reordering
		* ssno - state-space with no reordering
	* LP - all LP algorithms
	* DP - all dynamic programming algorithms
		* DPn - dynamic programming with no reodering
		* DPr - dynamic programming with random reodering
		* DPs - dynamic programming with 'sametogether' reodering
		* DPl - dynamic programming with 'leastout' reodering
		* DPnr - dynamic programming with no and random reoderings
		* DPsl - dynamic programming with 'sametogether' and 'leastout' reoderings
	* DPLPn - LP and dymaic programming with no reordering
* -r how many inputs to read from a single file (2 to 100)
* -min define lower bound of the range of numbers of wells to test
* -max define lower bound of the range of numbers of wells to test
