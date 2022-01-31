# pipette_opt

Latest change: 9 October 2021

This repository contains the Python 3.8 code implementing the algorithms for opitmising the pipette tip uptake during the distribution of DNA parts across the well array in which DNA assembly constructs are prepared. The detailed explanation of the LP algorithm can be found in the article 'A linear programming-based strategy to save pipette tips in automated DNA assembly' (K. Sechkar et al., 2021).

The code is organised as a distributable Python 3 package to be downloaded and used as a complement to DNA automation pipelines for Opentrons OT-2. The relevant python files contain both the API code and instructions for integrating our package with the supported automation pipelines (see below for details).  Currently, the [DNA-BOT](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) (BASIC assembly standard) and [OT2 Modular Cloning (MoClo) and Transformation in E.coli Workflow](https://github.com/DAMPLAB/OT2-MoClo-Transformation-Ecoli) (MoClo assembly standard) packages are supported.

The required Python packages and their versions are given in _requirements.txt_.

## Running test examples from the publication

In order to run the example considered in the publication, open the terminal, go to /src folder, and run the following command:
* python ppopt/test/testing.py -w LPnr -r 50 -min 2 -max 96 -n Start-Stop_Assembly_Random_Inputs.csv

This means that for every input size 2 &#8804; n &#8804; 96, first 50 inputs of this size listed in the file _Start-Stop_Assembly_Random_Inputs.csv_ are considered. 
The LP-based algorithm is run for each of these 50 inputs, with and without randomly permuting the read order of DNA parts;
then, the mean, the median and the standard deviation of the 50 optimised pipette tip consumptions are recorded.
Please refer to the 'Algorithm testing' section of this manual for a detailed explanation of how the program testing.py works and what its arguments are.

Completing the test run for all 95 input sizes takes >24h, so it is recommended to launch the program on a server,
rather than a laptop or a PC. 

To test the outcomes for the [DNA-BOT](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) package for automating BASIC assembly, download the DNA-BOT package (version 14 July 2020), and follow the instructions to modify the code.
Run 'dnabot_app.py', choosing 'storch_et_al_cons.csv' as the construct file and 'BIOLEGIO_BASIC_STD_SET_2_columns.csv' and 'part_plate_2_230419.csv' as part source files.

To test the outcomes for the [MoClo automation](https://github.com/DAMPLAB/OT2-MoClo-Transformation-Ecoli) package, download the DNA-BOT package (version 18 November 2020), and follow the instructions to modify the code.
Run 'moclo_transform_generator.py', choosing to make 'single' combinations, and selecting 'input-dna-map.csv' as the DNA plate map and 'combination-to-make-72.csv' as the list of DNA part combinations to assemble together. 

## Algorihm implementations
There are three subpackages, each of which implements one of the three approaches to solving the tip consumption optimisation problem:
* _lp_ - dividing the problem into a series of Linear Programming problems, and using the [Google OR-tools CP-SAT](https://developers.google.com/optimization/cp/cp_solver) or [GUROBI](https://www.gurobi.com/) optimiser to solve them. This algorithm is the one considered in the publication
* _statespace_ - searching a tree graph of states of the system (work in progress)
* _dp_ - dynamic programming (work in progress)

All of the algorithms receive an input in an abstract format independent of the assembly standard: each DNA part is represented as a tuple (a,b), i.e. species number b in the list of parts found on position a in the assembled contructs. The dictionary _caps_  outlines how many doses of each part's solution the pipette can hold; each nested array in the  2D-list _w_ outlines the composition of a single construct to be prepared. The output is an array _fin_, where each entry stands for the addition of a given DNA part to a single construct well; it also specifies if the tip must be changed to perform this operation.

The main() functions in these files allow to run the algorithms with certain inputs to test them. By changing and (un)commenting code lines in main, the algorithm can be run on a defined test input or on a randomly-generated input of up to 96 constructs.  Depending on what line of the code is uncommented, different algorithms from the same file can be tested out. By changing the reord argument, the reordering of the operation list (preprocessing of input to improve algorithm performance) can be selected.

## Adaptation for different DNA assemblies
The package enables pipette tip saving in existing automated DNA assembly protocols for Opentrons OT-2 for three different standards. The following files contain the API between the algorithms and the automation package, as well as instruction how to incorporate pipette tip saving into the automation package's workflow.
* _startstop.py_ contains API for the automated [Start-Stop](https://github.com/zoltuz/dna_assembler)  assembly standard.
* _basic.py_ contains API for the [DNA-BOT](https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT) package (BASIC assembly standard)
* _moclo.py_ contains APU for the package [OT2 Modular Cloning (MoClo) and Transformation in E.coli Workflow](https://github.com/DAMPLAB/OT2-MoClo-Transformation-Ecoli) package (MoClo assembly standard).

Please note that while the conversion of input into the internal abstract format and pipette tip consumption optimisation have been validated _in silico_, no real DNA assembly setup has been performed by a real Opentrons OT-2 robot in accordance with the obtained protocol.

For the Start-Stop and BASIC assembly, the resultant sequence of actions is saved in a .p file. By running the vis_cont.py programme and reading this .p file, the possible contaminations between consturct wells can be visualised.

## Algorithm testing
The subpackage _testing_ contains functions that allow to test the algorithms' performance.

The program _pytest.py_ uses the pytest function assert to validate the algorithms' performance on an input with all different parts (hence unoptimisable tip consumption) and a small 10-construct input with known minimum and maximum possible tip uptake.

Randomly generated inputs in the internal abstract format can be saved in .csv files:
* _input_generator.py_ contains the functions randomly generating Start-Stop assembly inputs with a given number of wells and parts of each type and saving them in .csv files
* _test_inp.py_ for every value of n from 2 to 96 generates 100 inputs, each describing n constructs;
  all these inputs are then saved in a csv file. 
  Change the main function to create test input collections for different ranges of n, different number of inputs per one n value,
  and different numbers of possible promoters/RBSs/CDSs/terminators.
  The file _Start-Stop_Assembly_Random_Inputs.csv_ was created using this function.

The program _testing.py_ reads the files created as test_inp.py and finds mean numbers of tips required and algorithm runtimes, along with standard deviations, saving them as .csv entries.
The program should be run from the command prompt with the following arguments
* -w: which algorithms in combination with which part reorderings to test:
	* all - all algorithms and all reorderings
	* LP - LP-based algorithm with all reordering 
	    * LPnr - LP-based algorithm with no or random part reodering
	* statespace - all state-space algorithms
		* sssametogether - state-space with 'sametogether' reordering
		* ssno - state-space with no reordering
	* DP - all dynamic programming algorithms
		* DPn - dynamic programming with no reodering
		* DPr - dynamic programming with random reodering
		* DPs - dynamic programming with 'sametogether' reodering
		* DPl - dynamic programming with 'leastout' reodering
		* DPnr - dynamic programming with no and random reoderings
		* DPsl - dynamic programming with 'sametogether' and 'leastout' reoderings
	* DPLPn - LP and dymaic programming with no reordering
* -r: how many inputs to read from a single file (2 to 100)
* -min: define lower bound of the range of numbers of wells to test
* -max: define lower bound of the range of numbers of wells to test

Note: the publication only considers the 'no reordering' and 'random part reordering' options.
None of the other reorderings have been tested and shown to consistently improve the algorithms' performance.

## Runtime estimation
By running _estimate_runtimes.py_, it is possible to use the simulation feature for the new Opentrons OT-2 API (version 2) to estimate the time that executing a program with Opentrons OT-2 would take. Please refer to the comments at the beginning of the script to select the value of WHICH_TEST that enables the simulation of the desired program.

Please note that this script is intended for runtime estimation only, and is not compatible with the rest of the package, which uses the old Opentrons OT-2 API version 1. Thus, opentrons 4.7.0, and not opentrons 3.21.0, is required to run this program.
