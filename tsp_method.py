# TSP-BASED METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.1.0, 22.7.20

# The project makes use of the 'tspy' package

import numpy as np
import time

# import functions from own files
from input_generator import wgenerator
from auxil import *
from tsp_reorder import leastout, sametogether, reorder_iddfs, reorder_greedy, reorder_a_star
from tsp_lp_solver import *

from tspy import TSP  # TSP solver package
from tspy.solvers.utils import get_cost
from tspy.solvers import TwoOpt_solver

# -------------------------------CLASS DEFINITIONS-------------------------------
# Each reagent matched with a subest of wells it is added to
class Ss:
    def __init__(self, reag, wellno):  # initialisation
        self.reag = reag
        self.wells = [wellno]

    def nuwell(self, wellno):  # record new well in the subset
        self.wells.append(wellno)

    def __str__(self):  # for printing the subset's reagent type and wells out
        strRep = self.reag + '|'
        for i in range(0, len(self.wells)):
            strRep = strRep + ' ' + str(self.wells[i])
        return strRep


# Final output format is an array of Operations - will be the same for ALL methods
class Oper:
    def __init__(self, reag, well):
        self.reag = reag
        self.well = well

    def __str__(self):  # for printing the subset's reagent type and wells out
        strRep = self.reag + ' -> w' + str(self.well)
        return strRep


# -------------------------------INPUT-------------------------------
# Will be replaced by a test example generator or manual input reading function

# this is the example given to me in the main pipette_opt file
w = [['p1', 'r2', 'c4', 't2'],
     ['p1', 'r2', 'c1', 't2'],
     ['p1', 'r3', 'c2', 't1'],
     ['p2', 'r3', 'c1', 't1']]


# -------------------------------MAIN-------------------------------
def main():
    fin = []  # final array where the operations are to be recorded

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to define the number of wells
    change 4 last arguments to define the size of p, r, c and t reagent sets"""
    w = wgenerator(96, 6, 6, 3, 4)

    cap = capac(pipcap=300, dose=40, airgap=10)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # the actual solver. Input empty file name to have w as input, empty w to use a json file as input
    tsp_method(w, fin, reord=None, filename=None, cap=cap)

    dispoper(fin)

    # TEST ONLY: to compare with previous method (no capacity)
    fin1=[]
    tsp_method(w,fin1,reord='',filename=None,cap=None)


    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')
    #print('The total number of pipette tips used is ' + str(tips)) #IRRELEVANT WITH LP SOLVER => COMMENTED
    print('The total number of pipette tips used is (independent calculation) ' + str(route_cost_with_w(fin, w, cap)))


# ---------------------SOLVER FUNCTION-------------------------------
# solves the problem, returns total cost
def tsp_method(w, fin, reord, filename, cap):
    subsets = []  # array of all subsets (class Ss)
    tips = 0  # counts the total number of tip changes

    # get the subsets
    if (filename == None):
        w_to_subsets(w, subsets)
    else:
        dic = jsonreader(filename, subsets=subsets, w=None)

    D = np.zeros((len(w), len(w)))  # initialise the matrix of distances, i.e. our graph of wells
    for i in range(0, len(D)):  # forbid going from one node to itself by setting a very high cost
        D[i][i] = 1000 * len(D)

    # print subsets and D (TEST ONLY)
    # disp(subsets, D)

    # reorder the subsets...
    if (reord == 'random'):  # ...randomly
        np.random.shuffle(subsets)
    elif (reord == 'random with time seed'):  # ...randomly using time as a seed
        np.random.RandomState(seed=round(time.time())).shuffle(subsets)
    elif (reord == 'leastout'):  # ...leastout
        leastout(subsets, len(w))
    elif (reord == 'sametogether'):  # ...sametogether
        sametogether(subsets, len(w))
    elif(reord!=None):  # (various state-space reorderings)
        origsubs = subsets.copy()
        subsets = []
        if (reord == 'nearest neighbour'):  # ...nearest neighbour algorithm (i.e. iddfs depth 1)
            reorder_iddfs(origsubs, subsets, D.copy(), 1)
        elif (reord == 'iddfs depth 2'):  # ...iddfs
            reorder_iddfs(origsubs, subsets, D.copy(), 2)
        elif (reord == 'greedy'):  # ...greedy tree search
            reorder_greedy(origsubs, subsets, D.copy(), 'countall')
        elif (reord == 'a*'):  # ...A*  tree search
            reorder_a_star(origsubs, subsets, D.copy(), 'countall')

    # print subsets and D (TEST ONLY)
    # disp(subsets, D)

    # implement the algorithm
    for i in range(0, len(subsets)):
        singlesub(subsets[i], D, fin, cap)
        print(str(i+1) + ' of ' + str(len(subsets)) + ' subsets processed')


# ---------------------SUBSET DISPLAY-------------------------------
def disp(subsets, D):
    for i in range(0, len(subsets)):
        print(subsets[i])
    print(D)


# -------------------------------SOLVE TSP FOR ONE SUBSET-------------------------------
def singlesub(subset, D, fin, cap):
    # PART 1: initial preparations
    # get length to avoid calling len too often
    sublen = len(subset.wells)

    # initialise the subset's matrix subD
    subD = np.zeros((sublen + 1,
                     sublen + 1))  # vertex 0, all edges to and from it being zero, allows to use cyclic TSP soluction for our PATH problem
    subD[0][0]=len(D)*1000
    # PART 2: select submatrix and update D as if problem for the subset is already solved
    for i_well in range(0, sublen):
        current_well = 0
        for j_D in range(0, len(D)):
            if (j_D == subset.wells[current_well]):
                subD[i_well + 1][current_well + 1] = D[subset.wells[i_well]][
                    j_D]  # select the edges within the subset into the submatrix
                if (current_well < sublen - 1):
                    current_well += 1
            else:
                D[subset.wells[i_well]][
                    j_D] = 1  # make the edge going from the subset into the rest of D equal to one (updating D)

    # PART 3: solve TSP for the subset
    if(cap!=None): #adjusting for capacity
        if (len(subD) == 2):
            routes = [[0, 1]]
        else:
            routes = lp_cap(subD, cap)

        # record in fin
        for route in routes:
            for i in range(0,len(route)):
                fin.append(Oper(subset.reag, subset.wells[route[i]-1]))

    else: # no adjustment for capacity
        if (len(subD) == 2):
            routes = [[0, 1]]
        else:
            tour = tsp_lp_gurobi(subD)
        # record in fin
        for i in range(0,len(tour)):
            fin.append(Oper(subset.reag, subset.wells[tour[i]-1]))

    # PART 5: return the adjusted number of pipette tip changes [IRRELEVANT WITH LP SOLVER => COMMENTED]
    #return tips

# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()

