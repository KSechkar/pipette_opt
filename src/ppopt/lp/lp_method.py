# LP-BASED METHOD (VEHICLE ROUTING PROBLEM) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v1.0.0, 23.2.21

"""
Having received the array of wells and parts w,
all additions to be performed are grouped into subsets according to the part added.
These subsets are reordered according to the method specified by the reord argument.

Then, for each subset a Linear Programming problem is defined and solved on a subgraph of the oriented
graph whose edge costs tell if a pipette tip is necessary between two wells (described by the distance matrix D).

The resultant sequence of operations is recorded in a class Oper list fin.
"""

from ppopt.auxil import *
from .lp_reorder import reorder_nns, reorder_greedy
from .lp_solver import *


# ---------------------SOLVER FUNCTION----------------------------------
# solves the problem
def lp_method(w, fin, reord, caps, maxtime):
    # PART 1: initial preparations

    # PART 1.1: get the subsets
    # (each 'subset' contains all additions of a given part; number of parts = number of subsets)
    subsets = []  # array of all subsets (class Ss)
    w_to_subsets(w, subsets)

    # PART 1.2: get the matrix of distances for the graph of wells
    D = np.zeros((len(w), len(w)))  # initialise
    for i in range(0, len(D)):  # forbid going from a node to itself by setting a very high cost
        D[i][i] = 1000 * len(D)

    # PART 1.3: set default maximum optimisation time if none specified
    if(maxtime==None):
        maxtime=1


    # PART 2: reorder the subsets
    if(reord!=None):
        if (reord == 'random'):  # ...randomly
            np.random.shuffle(subsets)
        elif (reord == 'random with time seed'):  # ...randomly using time as a seed
            np.random.RandomState(seed=round(time.time())).shuffle(subsets)
        elif (reord == 'leastout'):  # ...leastout
            leastout(subsets, w)
        elif (reord == 'sametogether'):  # ...sametogether
            sametogether(subsets, w)
        elif((reord[0:3]=='nns') or (reord=='greedy')):  # (various state-space reorderings)
            origsubs = subsets.copy()
            subsets = []
            if (reord == 'nns'):  # ...nearest neighbour algorithm (i.e. nns depth 1)
                reorder_nns(origsubs, subsets, D.copy(), 1, caps)
            elif (reord == 'nns depth 2'):  # ...nns depth 2
                reorder_nns(origsubs, subsets, D.copy(), 2, caps)
            elif (reord == 'greedy'):  # ...greedy tree search
                reorder_greedy(origsubs, subsets, D.copy(), 'countall',caps)

    # PART 3: implement the algorithm
    for i in range(0, len(subsets)):
        #print(str(i)+' of '+str(len(subsets))+' parts')  # uncomment if need to track the progress
        # call single-subset LP solver for each subset
        singlesub(subsets[i], D, fin, caps[subsets[i].part],maxtime)

    # PART 4: fix redundant tip changes
    fix_redundant(fin,w,caps)


# ------------------SOLVER FOR ONE SUBSET-------------------------------
# creates and solves an LP problem for
def singlesub(subset, D, fin, cap,maxtime):
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

    # PART 3: solve LP problem for the subset

    # 3a): capacitated problem
    if(cap!=None):
        # get the chain coverage
        if (len(subD) == 2):
            chains = [[1]]
        else:
            chains = lp_cap(subD, cap, maxtime)

        # record operations in fin
        for chain in chains:
            fin.append(Oper(subset.part, subset.wells[chain[0] - 1]))
            fin[-1].changed=True
            for i in range(1,len(chain)):
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


# -----------------CHECK & FIX REDUNDANT TIP CHANGES---------------------------
# If the same part is in many wells, the LP method can yield unnecessary tip changes.
# This program fixes them by calculating tip changes of the same operation list independently

def fix_redundant(fin,w,caps):
    # PART 1: initial preparations

    # make a copy, reset its tip change indicators
    cfin = deepcopy(fin)
    for i in range(0, len(cfin)):
        cfin[i].changed = False

    # create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])), dtype=bool)  # added[i][j]==True if part at address j in well i has been added

    # PART 2: get the cost

    # PART 2.1: the first operation in cfin
    cost = 1  # beginning the distribuiton => new tip taken
    cfin[0].changed = True  # indicate the tip's been changed
    added[cfin[0].well][cfin[0].part[0]] = 1  # indicate the part's been added

    # PART 2.2: all other operations
    for i in range(1, len(cfin)):
        one_cost = cost_func_with_w(cfin[0:i], cfin[i], w, added, caps)  # get operation cost
        cost += one_cost  # add operation cost

        added[cfin[i].well][cfin[i].part[0]] = 1  # indicate the part's been added
        # if the tip's been changed (operation cost 1), indicate that
        if (one_cost == 1):
            cfin[i].changed = True

    # PART 3: if the independently-determined tip changes are better, use this better operation list instead
    if(cost<route_cost(fin)):
        for i in range(0,len(fin)):
            fin[i]=cfin[i]

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
    w = [[(0, 1), (1, 2), (2, 4), (3, 1)],
         [(0, 2), (1, 2), (2, 1), (3, 1)],
         [(0, 1), (1, 2), (2, 2), (3, 2)],
         [(0, 2), (1, 3), (2, 1), (3, 1)]]

    w = wgenerator(96, 6, 6, 3, 4)

    # generate required volumes (for testing). Values taken from a real instance of Start-Stop assembly
    ss=[]
    w_to_subsets(w,ss)
    reqvols = {}
    for s in ss:
        if(s.part[0]==0):
            reqvols[s.part]=1.09
        elif(s.part[0]==1):
            reqvols[s.part]=0.33
        elif (s.part[0] == 2):
            reqvols[s.part] = 0.36
        else:
            reqvols[s.part] = 0.75

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # get capacitites
    caps=capacities(reqvols,10,1.0)

    # Call the solver. Specify the heuristic reordering used by changing reord; specify maximum optimisation time by changing maxtime
    lp_method(w, fin, reord='sametogether', caps=caps, maxtime=1)

    # display the solution
    dispoper(fin)

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')
    print('The total number of pipette tips used is (from resultant list) ' + str(route_cost(fin)))
    print('The total number of pipette tips used is (independent calculation) ' + str(independent_cost(fin, w, caps)[0]))

# main call
if __name__ == "__main__":
    # import necessary modules for independent test running
    import numpy as np
    import time
    from src.ppopt import wgenerator
    main()

