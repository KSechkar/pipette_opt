# CHECKING BASIC FUNCTIONS BY PYETEST
# By Kirill Sechkar
# v0.0.1, 11.9.20

from input_generator import wgenerator
from auxil import *
from lp_method import lp_method
from statespace_methods import nns, greedy_tree

# PART 1: check with a 3-part input of 96 wells wehere all parts are different
# The number of pipette tips used must be the maximum in all cases

# function creating the input
def make_all_different():
    w=[]
    typenames = 'prctabdefghijklmnoqsuvwxyz'
    for i in range(0, 96):
        well_entry = []
        for j in range(0, 3):
            well_entry.append(typenames[j] + str(i))
        w.append(well_entry.copy())
    return w

def check_all_diff():
    print('Checking the situation with all different parts...\n')

    # get input
    w_all_diff=make_all_different()

    # create capacity information (required volume 1 for all)
    ss = []
    w_to_subsets(w_all_diff, ss)
    reqvols = {}
    for s in ss:
        reqvols[s.part] = 1
    caps = capacities(reqvols, 10, 1.0)

    # for LP method
    fin=[]
    lp_method(w_all_diff,fin,None,caps)
    assert (independent_cost(fin,w_all_diff,caps) == (3*96))

    # for state-space methods
    # nns
    fin=[]
    nns(w_all_diff,fin,1,'sametogether',caps)
    assert (independent_cost(fin,w_all_diff,caps) == (3*96))

    # greedy search
    fin=[]
    greedy_tree(w_all_diff,fin,'optimistic+cap',None,caps)
    assert (independent_cost(fin,w_all_diff,caps) == (3*96))

# call the checker function
check_all_diff()

# PART 2: test with a known 4-part 10-well input
w_known = [['p0', 'r0', 'c0', 't0'],
           ['p0', 'r1', 'c1', 't1'],
           ['p0', 'r2', 'c2', 't2'],
           ['p0', 'r3', 'c3', 't3'],
           ['p0', 'r4', 'c4', 't4'],
           ['p0', 'r4', 'c5', 't5'],
           ['p0', 'r4', 'c6', 't6'],
           ['p1', 'r4', 'c7', 't7'],
           ['p2', 'r4', 'c8', 't8'],
           ['p3', 'r5', 'c9', 't9']]

def check_known():
    print('Checking the situation with a known input...\n')

    # get capacities
    ss = []
    w_to_subsets(w_known, ss)
    reqvols = {}
    for s in ss:
        reqvols[s.part] = 1
    caps = capacities(reqvols, 10, 1.0)

    # check if the number of subsets is right
    assert (len(ss) == 30)

    # check for every method
    methods = ['Nearest Neighbour', 'nns depth 2', 'Greedy',
               'Nearest Neighbour+sametogether', 'nns depth 2+sametogether', 'Greedy+sametogether',
               'LP+random', 'LP', 'LP+sametogether', 'LP+leastout',
               'LP+nearest neighbour', 'LP+nns depth 2', 'LP+greedy']
    for m in methods:
        # display the current method under investigation
        print(m)
        
        # solve
        fin = []
        if (m[0:2] == 'LP'):
            if (len(m) == 2):
                lp_method(w_known, fin, reord=None, caps=caps)
            else:
                lp_method(w_known, fin, m[3:], caps=caps)
        else:
            # define reordering
            if (m[0][-12:] == 'sametogether'):
                reord = 'sametogether'
            else:
                reord = None

            # get solution
            if (m[:7] == 'Nearest'):
                nns(w_known, fin, 1, reord, caps)
            elif (m[:3] == 'nns'):
                nns(w_known, fin, 2, reord, caps)
            elif (m[:6] == 'Greedy'):
                greedy_tree(w_known, fin, 'optimistic+cap', reord, caps)

        # the output must contain ALL operations
        assert(len(fin)==40)
        # get route cost and record
        cost = independent_cost(fin, w_known, caps)
        # the cost must be between the minimum and maximum possible values (which are known)
        assert ((cost <= 40) and (cost >= 31))
    
# call function
check_known()