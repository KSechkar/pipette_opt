# CHECKING BASIC FUNCTIONS BY PYETEST
# By Kirill Sechkar
# v0.1.0, 9.4.21

from auxil import *
from lp_method import lp_method
from dp_method import dp_method
from statespace_methods import nns, greedy_tree


# -----------------FUNCTION DEFINITIONS------------------
# function creating the input of wells (numwells in total) with all-different parts (numparts in each construct)
def all_different_input(numwells,numparts):
    w=[]
    for i in range(0, numwells):
        well_entry = []
        for j in range(0, numparts):
            well_entry.append((j,i))
        w.append(well_entry.copy())
    return w


# function solving the problem for well set w and capacity set caps with all methods specified in methods; 
# resultant costs, calculated from output sequence and independently, are returned as a dictionary
# (method name: [cost from fin, independent cost])
def validate_each(w, caps, methods):
    costs={}
    
    # solve and validate for each method
    for method in methods:
        print('Checking '+method+'...')
        fin = [] # empty array for output
        if (method[0:2] == 'LP'): # LP methods
            if (len(method) == 2):
                lp_method(w, fin, reord=None, caps=caps, maxtime=1)
            else:
                lp_method(w, fin, method[3:], caps=caps, maxtime=1)

        elif (method[0:2] == 'DP'): # dynamic programming methods
            if (len(method) == 2):
                dp_method(w, fin, reord=None, caps=caps)
            else:
                dp_method(w, fin, method[3:], caps=caps)

        else: # state-space methods
            # define reordering
            if (method[-12:] == 'sametogether'):
                reord = 'sametogether'
            else:
                reord = None
            # get solution
            if (method[:7] == 'Nearest'):
                nns(w, fin, 1, reord, caps)
            elif (method[:3] == 'NNs'):
                nns(w, fin, 2, reord, caps)
            elif (method[:6] == 'Greedy'):
                greedy_tree(w, fin, 'optimistic+cap', reord, caps)
        costs[method]=[route_cost(fin),independent_cost(fin, w, caps)[0]]
    
    # return the dictionary of costs
    return costs


# --------------------MAIN FUNCTION----------------------
def main():
    # array containing descriptors of all method tested
    methods = ['Nearest Neighbour', 'NNs depth2', 'Greedy',
               'Nearest Neighbour+sametogether', 'NNs depth2+sametogether', 'Greedy+sametogether',
               'Nearest Neighbour+leastout', 'NNs depth2+leastout', 'Greedy+leastout',
               'LP', 'LP+random', 'LP+sametogether', 'LP+greedy', 'LP+nearest neighbour', 'LP+nns depth 2',
               'DP', 'DP+random', 'DP+sametogether', 'DP+leastout']

    # PART 1: check with a 5-part input of 48 wells where all parts are different
    # The number of pipette tips used must be the maximum possible in all cases
    print('Checking the situation with all different parts...\n')

    # get input
    w_alldiff = all_different_input(48,5)

    # create capacity information (required volume 1, air gap 2 for all)
    ss_alldiff = []
    w_to_subsets(w_alldiff, ss_alldiff)
    reqvols = {}
    for s in ss_alldiff:
        reqvols[s.part] = 1
    caps_alldiff = capacities(reqvols, 10, 2.0)
    
    # check if the number of subsets is right
    assert (len(ss_alldiff) == 48*5)

    # check if all capacities are indeed 4 like we would expect
    # (maximum arrangement is: solution-gap-solution-gap-solution-gap-solution, 10 uL in total)
    for cap in caps_alldiff.values():
        assert (cap == 4)

    # run the algorithm for all methods
    costs_alldiff=validate_each(w=w_alldiff, caps=caps_alldiff, methods=methods)
    
    #  compare to the expected correct cost of maximum value 48*5
    for costpair in costs_alldiff.values():
        assert (costpair[0] == 48 * 5)
        assert (costpair[1] == 48 * 5)

    # PART 2: test with a known 3-part 10-well input
    print('Checking the situation with a known input...\n')

    # known input definition
    w_known = [[(0,0), (1,0), (2,0), (3,0)],
               [(0,0), (1,1), (2,1), (3,1)],
               [(0,0), (1,2), (2,2), (3,2)],
               [(0,0), (1,3), (2,3), (3,3)],
               [(0,0), (1,4), (2,4), (3,4)],
               [(0,0), (1,4), (2,5), (3,5)],
               [(0,0), (1,4), (2,6), (3,6)],
               [(0,1), (1,4), (2,7), (3,7)],
               [(0,2), (1,4), (2,8), (3,8)],
               [(0,3), (1,5), (2,9), (3,9)]]

    # get capacities (required volume 1.5, air gap 1 for all)
    ss_known = []
    w_to_subsets(w_known, ss_known)
    reqvols = {}
    for s in ss_known:
        reqvols[s.part] = 1.5
    caps_known = capacities(reqvols, 10, 1.0)

    # check if the number of subsets is right
    assert (len(ss_known) == 30)

    # check if all capacities are indeed 4 like we would expect
    # (maximum arrangement is: solution-gap-solution-gap-solution-gap-solution, 9 uL in total)
    for cap in caps_known.values():
        assert (cap == 4)

    # run the algorithm for all methods
    costs_known = validate_each(w=w_known, caps=caps_known, methods=methods)

    # the cost must be between the minimum and maximum possible values (which are known)
    for costpair in costs_known.values():
        assert ((costpair[0] >= 31) and (costpair[0] <= 40))
        assert ((costpair[0] >= 31) and (costpair[0] <= 40))


# main call
if __name__ == "__main__":
    main()