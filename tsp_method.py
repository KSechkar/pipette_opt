# TSP-BASED METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.1.0, 22.7.20

"""
Having received the array of wells and parts w,
all additions to be performed are grouped into subsets according to the part added.
These subsets are reordered according to the method specified by the reord argument.

Then, for each subset a Linear Programming problem is defined and solved on a subgraph of the oriented
graph whose edge costs tell if a pipette tip is necessary between two wells (described by the distance matrix D).

The resultant sequence of operations is recorded in a class Oper list fin.
"""

import numpy as np
import time

# import functions from own files
from input_generator import wgenerator
from auxil import *
from tsp_reorder import leastout, sametogether, reorder_nns, reorder_greedy
from tsp_lp_solver import *

# ----------------------CLASS DEFINITIONS-------------------------------
# Each part matched with the wells it is added to
class Ss:
    # initialisation
    def __init__(self, part, wellno):
        self.part = part
        self.wells = [wellno]

    # record new well in the subset
    def nuwell(self, wellno):
        self.wells.append(wellno)

    # for printing the subset's part type and wells
    def __str__(self):
        strRep = self.part + '|'
        for i in range(0, len(self.wells)):
            strRep = strRep + ' ' + str(self.wells[i])
        return strRep


# Final output format is an array of Operations - the same for ALL methods
class Oper:
    # initialisation
    def __init__(self, part, well):
        self.part = part
        self.well = well
        self.changed = False

    # for printing the part type and destination well
    def __str__(self):
        strRep = self.part + ' -> w' + str(self.well)
        return strRep


# ---------------------SOLVER FUNCTION----------------------------------
# solves the problem
def tsp_method(w, fin, reord, caps):
    # PART 1: initial preparations

    # PART 1.1: get the subsets
    # (each 'subset' contains all additions of a given part; number of parts = number of subsets)
    subsets = []  # array of all subsets (class Ss)
    w_to_subsets(w, subsets)

    # PART 1.2: get the matrix of distances for the graph of wells
    D = np.zeros((len(w), len(w)))  # initialise
    for i in range(0, len(D)):  # forbid going from a node to itself by setting a very high cost
        D[i][i] = 1000 * len(D)


    # PART 2: reorder the subsets
    if (reord == 'random'):  # ...randomly
        np.random.shuffle(subsets)
    elif (reord == 'random with time seed'):  # ...randomly using time as a seed
        np.random.RandomState(seed=round(time.time())).shuffle(subsets)
    elif (reord == 'leastout'):  # ...leastout
        leastout(subsets, w)
    elif (reord == 'sametogether'):  # ...sametogether
        sametogether(subsets, w)
    elif(reord!=None):  # (various state-space reorderings)
        origsubs = subsets.copy()
        subsets = []
        if (reord == 'nearest neighbour'):  # ...nearest neighbour algorithm (i.e. nns depth 1)
            reorder_nns(origsubs, subsets, D.copy(), 1, caps)
        elif (reord == 'nns depth 2'):  # ...nns depth 2
            reorder_nns(origsubs, subsets, D.copy(), 2, caps)
        elif (reord == 'greedy'):  # ...greedy tree search
            reorder_greedy(origsubs, subsets, D.copy(), 'countall',caps)

    # PART 3: implement the algorithm
    for i in range(0, len(subsets)):
        # call single-subset LP solver for each subset
        singlesub(subsets[i], D, fin, caps[subsets[i].part])


# ------------------SOLVER FOR ONE SUBSET-------------------------------
# creates and solves an LP problem for
def singlesub(subset, D, fin, cap):
    # PART 1: initial preparations
    # get length to avoid calling len too often
    sublen = len(subset.wells)

    # initialise subD, the distance matrix for the subgraph
    subD = np.zeros((sublen + 1,sublen + 1)) # an extra 0 node is needed for tspy compatibility
    subD[0][0]=len(D)*1000

    # PART 2: select the submatrix and update D
    """
    For every well belonging to the subset (given by subset.wells[i_well]), we consider all outgoing edges
    (given by D[subset.wells[i_well]][j_D]).
    If the edge arrives at another subset well (given by subset.wells[current_well]), we select it into the subgraph,
    i.e. we get its value into subD[i_well+1][current_well].
    If the edge arrives into a well not in the subset, its cost in D must be updated to be 1.
    """
    for i_well in range(0, sublen):
        current_well = 0
        for j_D in range(0, len(D)):
            if (j_D == subset.wells[current_well]):
                subD[i_well + 1][current_well + 1] = D[subset.wells[i_well]][j_D] # selecting the edge
                if (current_well < sublen - 1):
                    current_well += 1
            else:
                D[subset.wells[i_well]][j_D] = 1  # updating D

    # PART 3: solve TSP for the subset

    # 3a): capacitated problem
    if(cap!=None):
        # get the chain coverage
        if (len(subD) == 2):
            chains = [[1]]
        else:
            chains = lp_cap(subD, cap,maxtime=None)

        # record operations in fin
        for chain in chains:
            for i in range(0,len(chain)):
                fin.append(Oper(subset.part, subset.wells[chain[i]-1]))

    # 3b): non-capacitated problem
    else:
        # get the TSP tour
        if (len(subD) == 2):
            tour = [0,1]
        else:
            tour = tsp_lp_gurobi(subD)

        # record operations in fin
        for i in range(1,len(tour)):
            fin.append(Oper(subset.part, subset.wells[tour[i]-1]))


# ---------------------DISPLAYING SUBSETS-------------------------------
def disp(subsets, D):
    for i in range(0, len(subsets)):
        print(subsets[i])
    print(D)

# -----------------------MAIN (TESTING ONLY)----------0-----------------
def main():
    fin = []  # final array where the operations are to be recorded

    """
    INPUT:
    a) Use a manually defined well array, or
    b) Generate 
            change 1st argument of wgenerator to define the number of wells/constructs
            change 4 last arguments of wgenerator to define the size of p, r, c and t part sets
    Comment out the respective line to deselect
    """
    w = [['p1', 'r2', 'c4', 't1'],
         ['p2', 'r2', 'c1', 't1'],
         ['p1', 'r3', 'c2', 't2'],
         ['p2', 'r3', 'c1', 't1']]

    w = wgenerator(96, 6, 6, 3, 4)

    # generate required volumes (for testing)
    ss=[]
    w_to_subsets(w,ss)
    reqvols = {}
    for s in ss:
        if(s.part[0]=='p'):
            reqvols[s.part]=1.09
        elif(s.part[0]=='r'):
            reqvols[s.part]=0.33
        elif (s.part[0] == 'c'):
            reqvols[s.part] = 0.36
        else:
            reqvols[s.part] = 0.75

    # get capacitites
    caps=capacities(reqvols,10,1.0)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # Call the solver. Input empty file name to have w as input, empty w to use a json file as input
    tsp_method(w, fin, reord=None, caps=caps)

    dispoper(fin)

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')
    print('The total number of pipette tips used is (independent calculation) ' + str(route_cost_with_w(fin, w, caps)))


# main call
if __name__ == "__main__":
    main()

