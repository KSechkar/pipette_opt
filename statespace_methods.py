# STATE SPACE-BASED METHODS (SHORTEST PATH/TREE SEARCH) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.1.0, 22.7.20

import numpy as np
import time

# import functions from own files
from input_generator import wgenerator
from auxil import *
from tsp_reorder import sametogether

# match reagents with their addresses in w
ADDRESS = {'p': 0, 'r': 1, 'c': 2, 't': 3}


# -------------------------------CLASS DEFINITIONS-------------------------------
# Final output format is an array of Operations - will be the same for ALL methods
class Oper:
    def __init__(self, reag, well):
        self.reag = reag
        self.well = well

    def __str__(self):  # for printing the subset's reagent type and wells out
        strRep = self.reag + ' -> w' + str(self.well)
        return strRep


# needed to construct a state space
class State:
    def __init__(self, added, trace, todo, g, f):
        self.todo = todo
        self.trace = trace
        self.added = added
        self.g = g
        self.f = f

    def __str_self(self):
        strRep = 'Last performed: ' + str(self.last) + '\n' + str(self.added)
        return strRep


# -------------------------------INPUT-------------------------------
# Will be replaced by a test example generator or manual input reading function

# this is the example given to me in the main pipette_opt file
w = [['p1', 'r2', 'c4', 't2'],
     ['p2', 'r2', 'c1', 't2'],
     ['p1', 'r3', 'c2', 't1'],
     ['p2', 'r3', 'c1', 't1']]


# -------------------------------MAIN-------------------------------
def main():
    fin = []  # an Oper list of operations in the order which they should be performed

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to dfine the number of wells
    change 4 last arguments to define the size of p, r, c and t reagent sets"""
    w = wgenerator(96, 6, 6, 3, 4)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # determine how many vector doses a pipette can hold (capacity)
    # for now, values selected manually in the order: pipette capacity, reagent volume, air gap volume
    cap = capac(pipcap=300,dose=40,airgap=10)

    # use nearest-neighbour tree search to solve the problem
    # iddfs(w, fin, 1,reord=None, cap=cap)

    # use iddfs to solve the problem
    iddfs(w, fin, 2, reord=None, cap=cap)

    # use a_star on a tree to solve the problem (NOT WORKING)
    # a_star_tree(w,fin,'optimistic')

    # use greedy algorithm on a tree to solve the problem
    # greedy_tree(w, fin, 'optimistic+cap', reord=None, cap=cap)

    dispoper(fin)

    # PERFORMANCE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')

    print('The total number of pipette tips used is ' + str(route_cost_with_w(fin, w,cap)))


# -------------------------------SOLVERS (IDDFS)-------------------------------
# iddfs function
def iddfs(w, fin, depth, reord,cap):
    ops = []  # an Oper list of operations to be performed
    getops(w, ops, reord)

    # if we want to randomise the operation order
    # np.random.shuffle(ops)

    all_operations = len(w) * len(w[0])

    fin.append(ops[0])
    ops.pop(0)

    # if (with_w):
    added = np.zeros((len(w), len(w[0])))  # tells which reagents were added to which well
    added[fin[0].well][ADDRESS[fin[0].reag[0]]] = 1

    while (len(fin) < all_operations):
        # if (with_w):
        nextop = iddfs_oneiter_with_w(ops, fin, 1, depth, w, added,cap)
        added[ops[nextop].well][ADDRESS[ops[nextop].reag[0]]] = 1
        # else:
            # nextop = iddfs_oneiter(ops, fin, 1, depth)
        fin.append(ops[nextop])
        ops.pop(nextop)


"""
# single iteration of iddfs
def iddfs_oneiter(ops, fin, curdepth, depth):
    # determine the potential cost of each possible operation
    potcost = []
    for i in range(0, len(ops)):
        potcost.append(cost_func(fin, ops[i]))
        # next iteration
        if (curdepth < depth and len(ops) != 1):
            potcost[i] += iddfs_oneiter(ops[0:i] + ops[(i + 1):len(ops)], fin + [ops[i]], curdepth + 1, depth)

    # act according to the determined costs
    if (curdepth == 1):
        answer = potcost.index(min(potcost))
    else:
        answer = min(potcost)
    return answer
"""


# single iteration of iddfs (with w)
def iddfs_oneiter_with_w(ops, fin, curdepth, depth, w, added,cap):
    # determine the potential cost of each possible operation
    potcost = []
    for i in range(0, len(ops)):
        potcost.append(cost_func_with_w(fin, ops[i], w, added,cap))
        # next iteration
        if (curdepth < depth and len(ops) != 1):
            # change the inputs for next iteration
            added[ops[i].well][ADDRESS[ops[i].reag[0]]] = 1
            fin.append(ops[i])
            ops.pop(i)

            # call next iteration
            potcost[i] += iddfs_oneiter_with_w(ops, fin, curdepth + 1, depth, w, added, cap)

            # change the inputs back
            ops.insert(i, fin[len(fin) - 1])
            fin.pop()
            added[ops[i].well][ADDRESS[ops[i].reag[0]]] = 0

    # act according to the determined costs
    if (curdepth == 1):
        answer = potcost.index(min(potcost))
    else:
        answer = min(potcost)

    return answer

"""
# -------------------------------SOLVER (A*)-------------------------------
# A* solver - on a tree
def a_star_tree(w, fin, heur, reord):
    ops = []  # an Oper list of operations to be performed
    getops(w, ops, reord)
    alloperations = len(ops)

    # if we want to randomise the operation order
    # np.random.shuffle(ops)

    # starting position - state with no operations performed
    states = [[]]  # list of states in state-space under consideration
    unstates = [ops]  # list of reagents NOT added for a given state
    g = [0]  # distance from origin, mirrors states
    f = [h_tree([], ops, heur)]  # f(states[i])=g(states[i])+h(states[i]), mirrors states
    l = [0]

    while (len(states) != 0):
        consider = f.index(min(f))

        # TEST ONLY
        #print(str(len(states)) + ' ' + str(h_tree(states[consider], unstates[consider], 'optimistic')) + ' ' + str(max(l)))

        # print(consider)
        if (len(states[consider]) == alloperations):
            for s in states[consider]:
                fin.append(s)
            break

        # add neighbours to considered states
        for i in range(0, len(unstates[consider])):
            states.append(states[consider] + [unstates[consider][i]])
            unstates.append(unstates[consider][0:i] + unstates[consider][i + 1:len(unstates[consider])])
            g.append(g[consider] + cost_func(states[consider], unstates[consider][i]))
            f.append(g[-1] + h_tree(states[-1], unstates[-1], heur))
            l.append(len(states[-1]))  # TEST ONLY

        # remove the state in question
        states.pop(consider)
        unstates.pop(consider)
        g.pop(consider)
        f.pop(consider)
        l.pop(consider)  # TEST ONLY
"""


# ---------------------------SOLVER (GREEDY)-------------------
def greedy_tree(w, fin, heur, reord,cap):
    ops = []  # an Oper list of operations to be performed
    getops(w, ops, reord)

    # if we want to randomise the operation order
    # np.random.shuffle(ops)

    all_operations = len(w) * len(w[0])

    fin.append(ops[0])
    ops.pop(0)

    added = np.zeros((len(w), len(w[0])))  # tells which reagents were added to which well
    added[fin[0].well][ADDRESS[fin[0].reag[0]]] = 1

    while (len(fin) < all_operations):
        #print(len(fin))
        nextop = greedy_tree_onestep(ops, fin, w, added, heur,cap)
        added[ops[nextop].well][ADDRESS[ops[nextop].reag[0]]] = 1

        fin.append(ops[nextop])
        ops.pop(nextop)


def greedy_tree_onestep(ops, fin, w, added, heur,cap):
    potcost = []
    for i in range(0, len(ops)):
        potcost.append(cost_func_with_w(fin, ops[i], w, added,cap))  # cost function
        fin.append(ops[i])
        ops.pop(i)
        potcost[-1] += h_tree(fin, ops, heur,cap)  # heurstic
        ops.insert(i, fin[-1])
        fin.pop()

    # act according to the determined costs
    return potcost.index(min(potcost))


# heuristic function
def h_tree(state, unstate, heur,cap):
    if (heur == 'optimistic'):  # h is the forecoming number of tip changes assuming we don't have to do extra changes
        already = []
        for unop in unstate:
            ispresent = False
            for alread in already:
                if (unop.reag == alread):
                    ispresent = True
                    break
            if not ispresent:
                already.append(unop.reag)
        return len(already)
    elif (heur == 'optimistic+cap'):
        subsets=[]
        ops_to_subsets(unstate,subsets)
        est_cost=0
        for subset in subsets:
            est_cost+=1
            est_cost+=round(len(subset.wells)/cap)
        return est_cost


# -------------------------------AUXILIARY FUNCTIONS-------------------------------
# get a list of all operations from w
def getops(w, ops, reord):
    # as iddfs and greedy search pick the FIRST element with minimum cost, the pre-set order matters, hence 'reord'
    if (reord == None):
        for well in range(0, len(w)):
            for reagent in range(0, len(w[well])):
                ops.append(Oper(w[well][reagent], well))
        return
    elif (reord == 'random'):
        for well in range(0, len(w)):
            for reagent in range(0, len(w[well])):
                ops.append(Oper(w[well][reagent], well))
        np.random.shuffle(ops)
        return
    if (reord == 'justsubsets' or reord == 'sametogether'):
        subsets = []
        w_to_subsets(w, subsets)
        if (reord == 'sametogether'):
            sametogether(subsets, len(w))
        subsets_to_ops(subsets, ops)


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
