# pipette_opt

This branch introduces heuristic reorderings (instead pf simply random) to our method

LEASTOUT
First, we go through the subsets with the least number of edges going from the subset to the remaining wells. Therefore, in the beginning we make the least number of edges equal to 1

SAMETOGETHER
Only one reagent of type P, R, T or C is added to a well. Thus, if we add all reagents of same type consequently, there will be no 'interaction' between the subsets (no 'oned' edges created by the same subset). The four blocks are then rearranged as in LEASTOUT
