# pipette_opt
Latest change: 4 May 2020 - improved comments, raw version of SAMETOGETHER added

This branch introduces a counter of pipette tip changes and heuristic reorderings (instead of simply random) to our method

COUNTING TIP CHANGES
singlesub() now returns the updated number of tip changes

LEASTOUT REORDERING
First, we go through the subsets with the least number of edges going from the subset to the remaining wells. Therefore, in the beginning we make the least number of edges equal to 1

SAMETOGETHER REORDERING (written, but may need additional testing)
Only one reagent of type P, R, T or C is added to a well. Thus, if we add all reagents of same type consequently, there will be no 'interaction' between the subsets (no 'oned' edges created by the same subset). The four blocks are then rearranged as in LEASTOUT
