# TSP-BASED METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.6.1, 10.7.20

# The project makes use of the 'tspy' package

import numpy as np
import time
from tspy import TSP  # TSP solver package
from tspy.solvers.utils import get_cost
from tspy.solvers import TwoOpt_solver

# import functions from own files
from input_generator import wgenerator
from auxil import dispoper, route_cost_with_w, jsonreader,w_to_subsets
from tsp_reorder import leastout, sametogether, reorder_iddfs, reorder_greedy, reorder_a_star


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
     ['p2', 'r2', 'c1', 't2'],
     ['p1', 'r3', 'c2', 't1'],
     ['p2', 'r3', 'c1', 't1']]


# w=[['p3', 'r5', 'c2', 't0'], ['p5', 'r2', 'c0', 't1'], ['p1', 'r5', 'c2', 't0'], ['p1', 'r5', 'c2', 't1'], ['p3', 'r3', 'c2', 't1'], ['p2', 'r5', 'c0', 't3'], ['p1', 'r4', 'c2', 't2'], ['p0', 'r0', 'c1', 't3'], ['p0', 'r1', 'c1', 't2'], ['p3', 'r4', 'c2', 't3'], ['p1', 'r0', 'c1', 't1'], ['p5', 'r3', 'c2', 't1'], ['p4', 'r0', 'c1', 't1'], ['p5', 'r0', 'c0', 't1'], ['p4', 'r5', 'c0', 't1'], ['p1', 'r0', 'c1', 't3'], ['p5', 'r4', 'c1', 't0'], ['p0', 'r3', 'c0', 't3'], ['p1', 'r1', 'c2', 't1'], ['p5', 'r3', 'c2', 't1'], ['p5', 'r2', 'c1', 't1'], ['p4', 'r3', 'c0', 't1'], ['p3', 'r5', 'c0', 't1'], ['p2', 'r0', 'c1', 't1'], ['p5', 'r2', 'c0', 't2'], ['p3', 'r4', 'c2', 't1'], ['p1', 'r3', 'c1', 't0'], ['p5', 'r3', 'c0', 't3'], ['p0', 'r0', 'c0', 't3'], ['p2', 'r3', 'c1', 't3'], ['p2', 'r5', 'c1', 't3'], ['p3', 'r3', 'c0', 't1'], ['p5', 'r5', 'c1', 't3'], ['p3', 'r1', 'c0', 't1'], ['p3', 'r1', 'c1', 't1'], ['p4', 'r5', 'c2', 't2'], ['p0', 'r3', 'c1', 't2'], ['p1', 'r5', 'c1', 't1'], ['p0', 'r5', 'c1', 't0'], ['p2', 'r1', 'c2', 't3'], ['p5', 'r2', 'c1', 't2'], ['p2', 'r1', 'c1', 't2'], ['p0', 'r0', 'c0', 't0'], ['p5', 'r4', 'c1', 't2'], ['p4', 'r0', 'c0', 't0'], ['p2', 'r3', 'c2', 't2'], ['p3', 'r1', 'c0', 't1'], ['p5', 'r1', 'c2', 't0'], ['p1', 'r5', 'c1', 't0'], ['p0', 'r5', 'c1', 't3'], ['p3', 'r0', 'c2', 't2'], ['p2', 'r2', 'c0', 't0'], ['p0', 'r5', 'c2', 't0'], ['p1', 'r4', 'c0', 't3'], ['p2', 'r2', 'c0', 't2'], ['p1', 'r2', 'c2', 't0'], ['p0', 'r5', 'c1', 't0'], ['p3', 'r3', 'c2', 't0'], ['p2', 'r0', 'c0', 't2'], ['p2', 'r2', 'c0', 't2'], ['p0', 'r5', 'c1', 't3'], ['p3', 'r2', 'c2', 't0'], ['p4', 'r0', 'c1', 't1'], ['p2', 'r1', 'c0', 't0'], ['p3', 'r1', 'c1', 't2'], ['p4', 'r1', 'c1', 't2'], ['p4', 'r1', 'c1', 't1'], ['p1', 'r3', 'c1', 't0'], ['p4', 'r1', 'c2', 't3'], ['p0', 'r2', 'c1', 't3'], ['p2', 'r3', 'c2', 't3'], ['p5', 'r3', 'c1', 't3'], ['p4', 'r0', 'c1', 't0'], ['p4', 'r5', 'c1', 't1'], ['p0', 'r2', 'c1', 't0'], ['p0', 'r1', 'c1', 't0'], ['p5', 'r0', 'c0', 't3'], ['p5', 'r4', 'c0', 't0'], ['p5', 'r5', 'c2', 't2'], ['p2', 'r2', 'c1', 't2'], ['p1', 'r1', 'c1', 't1'], ['p2', 'r5', 'c1', 't1'], ['p1', 'r2', 'c0', 't1'], ['p5', 'r5', 'c0', 't1'], ['p2', 'r3', 'c0', 't1'], ['p4', 'r4', 'c0', 't2'], ['p5', 'r4', 'c1', 't1'], ['p0', 'r4', 'c0', 't3'], ['p3', 'r1', 'c1', 't1'], ['p4', 'r1', 'c2', 't2'], ['p5', 'r4', 'c0', 't0'], ['p0', 'r3', 'c1', 't3'], ['p0', 'r4', 'c2', 't1'], ['p2', 'r4', 'c2', 't2'], ['p2', 'r5', 'c2', 't1'], ['p3', 'r2', 'c2', 't3']]

# -------------------------------MAIN-------------------------------
def main():
    fin = []  # final array where the operations are to be recorded

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to define the number of wells
    change 4 last arguments to define the size of p, r, c and t reagent sets"""
    w = wgenerator(96, 6, 6, 3, 4)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # the actual solver. Input empty file name to have w as input, empty w to use a json file as input
    tips = tsp_method(w, fin, 'sametogether',filename=None)

    dispoper(fin)

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')
    print('The total number of pipette tips used is ' + str(tips))
    print('The total number of pipette tips used is (independent calculation)' + str(route_cost_with_w(fin, w)))


# ---------------------SOLVER FUNCTION-------------------------------
# solves the problem, returns total cost
def tsp_method(w, fin, reord,filename):
    subsets = []  # array of all subsets (class Ss)
    tips = 0  # counts the total number of tip changes

    # get the subsets
    if (filename == None):
        w_to_subsets(w, subsets)
    else:
        dic=jsonreader(filename, subsets)

    D = np.zeros((len(w), len(w)))  # initialise the matrix of distances, i.e. our graph of wells
    for i in range(0, len(D)):  # forbid going from one node to itself by setting a very high cost
        D[i][i] = 1000 * len(D)

    tips = len(subsets)  # anyhow, we have to change the tip between the different reagents and we have a tip at first

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
        #print(str(i) + ' of ' + str(len(subsets) - 1) + ' subsets processed')
        tips = singlesub(subsets[i], D, fin, tips)

    return tips


# ---------------------SUBSET DISPLAY-------------------------------
def disp(subsets, D):
    for i in range(0, len(subsets)):
        print(subsets[i])
    print(D)


# -------------------------------SOLVE TSP FOR ONE SUBSET-------------------------------
def singlesub(subset, D, fin, tips):
    # PART 1: initial preparations
    # get length to avoid calling len too often
    sublen = len(subset.wells)

    # initialise the subset's matrix subD
    subD = np.zeros((sublen + 1,
                     sublen + 1))  # vertex 0, all edges to and from it being zero, allows to use cyclic TSP soluction for our PATH problem

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
    tsp = TSP()
    tsp.read_mat(subD)

    two_opt = TwoOpt_solver(initial_tour='NN', iter_num=100)
    tour = tsp.get_approx_solution(two_opt)

    # PART 4: record the operations into the final output, 'unwrapping' the cycle arround the added zero node to create a path
    # find the position of the zero node in the tour
    i = 0
    while (tour[i] != 0):
        i += 1
    # record the part after the zero node
    i += 1
    while (i < len(tour) - 1):
        fin.append(Oper(subset.reag, subset.wells[tour[i] - 1]))
        i += 1
    # record the part 'before'the zero node
    i = 0
    while (tour[i] != 0):
        fin.append(Oper(subset.reag, subset.wells[tour[i] - 1]))
        i += 1

    # PART 5: return the adjusted number of pipette tip changes
    return tips + get_cost(tour, tsp)  # include the tour cost in the number of tip changes


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
