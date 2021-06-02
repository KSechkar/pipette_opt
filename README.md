# pipette_opt

Latest change: 8 April 2021

This repository contains the Python 3.8 code implementing the algorithms for opitmising the pipette tip uptake during the distribution of DNA parts across the well array in which DNA assembly constructs are prepared. The detailed explanation of the algorithms can be found in the article 'Algorithms for pipette tip saving in automated DNA assembly'.

The code is organised as a distributable Python 3 package.

## Algotihm implementations
There are three subpacjages, each of which implements one of the three approaches to solving the tip consumption optimisation problem:
* _statespace_ - searching a tree graph of states of the system
* _lp_ - dividing the problem into a series of Linear Programming problems, and using the [GUROBI](https://www.gurobi.com/) optimiser to solve them
* _dp_ - dynamic programming

All of the algorithms receive an input in an abstract format independent of the assembly standard: each DNA part is represented as a tuple (a,b), i.e. species number b in the list of parts found on position a in the assembled contructs. The dictionary _caps_  outlines how many doses of each part's solution the pipette can hold; each nested array in the  2D-list _w_ outlines the composition of a single construct to be prepared. The output is an array _fin_, where each entry stands for the addition of a given DNA part to a single construct well; it also specifies if the tip must be changed to perform this operation.

The main() functions in these files allow to run the algorithms with certain inputs to test them. By changing and (un)commenting code lines in main, the algorithm can be run on a defined test input or on a randomly-generated input of up to 96 constructs.  Depending on what line of the code is uncommented, different algorithms from the same file can be tested out. By changing the reord argument, the reordering of the operation list (preprocessing of input to improve algorithm performance) can be selected.

## Adaptation for different DNA assemblies
The package enables pipette tip saving in existing automated DNA assembly protocols for Opentrons OT-2 for three different standards. The following files contain the API between the algorithms and the automation package, as well as instruction how to incorporate pipette tip saving into the automation package's workflow.
* _startstop.py_ contains API for the automated [Start-Stop](https://github.com/zoltuz/dna_assembler)  assembly standard.
* _basic.py_ contains API for the [DNA-BOT](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) package (BASIC assembly standard)
* _moclo.py_ contains APU for the package [OT2 Modular Cloning (MoClo) and Transformation in E.coli Workflow](https://github.com/DAMPLAB/OT2-MoClo-Transformation-Ecoli) package (MoClo assembly standard).

Please note that while the conversion of input into the internal abstract format and pipette tip consumption optimisation have been validated, no DNA assemblyhas been executed and tested in accordance with the obtained protocol.

For the Start-Stop and BASIC assembly, the resultant sequence of actions is saved in a .p file. By running the vis_cont.py programme and reading this .p file, the possible contaminations between consturct wells can be visualised.

## Algorithm testing
The subpackage _testing_ contains functions that allow to test the algorithms' performance.

The program _pytest.py_ uses the pytest function assert to validate the algorithms' performance on an input with all different parts (hence unoptimisable tip consumption) and a small 10-construct input with known minimum and maximum possible tip uptake.

Randomly generated inputs in the internal abstract format can be saved in .csv files:
* _input_generator.py_ contains the functions randomly generating Start-Stop assembly inputs with a given number of wells and parts of each type and saving them in .csv files
* _test_inp.py_ for n from 2 to 96 saves a .csv file that contains 100 inputs, n constructs each

The program _testing.py_ reads the files created as test_inp.py and finds mean numbers of tips required and algorithm runtimes, along with standard deviations, saving them as .csv entries. The program should be run from command prompt with the following arguments
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
