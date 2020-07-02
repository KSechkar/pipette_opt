# STATE SPACE-BASED METHODS (SHORTEST PATH/TREE SEARCH) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.2, 2.7.20

import numpy as np
import time

# import functions from own files
from input_generator import wgenerator
from auxil import route_cost, cost_func, dispoper, cost_func_with_w, route_cost_with_w

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

    # use iddfs to solve the problem
    iddfs(w, fin, 2, True)

    dispoper(fin)

    # PERFORMANCE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')

    print('The total number of pipette tips used is ' + str(route_cost_with_w(fin, w)))


# -------------------------------SOLVERS-------------------------------
# iddfs function
def iddfs(w, fin, depth, with_w):
    ops = []  # an Oper list of operations to be performed
    getops(w, ops)

    # if we want to randomise the operation order
    np.random.shuffle(ops)

    all_operations = len(w) * len(w[0])

    fin.append(ops[0])
    ops.pop(0)

    if (with_w):
        added = np.zeros((len(w), len(w[0])))  # tells which reagents were added to which well
        added[fin[0].well][ADDRESS[fin[0].reag[0]]] = 1

    while (len(fin) < all_operations):
        print(len(fin))
        if (with_w):
            nextop = iddfs_oneiter_with_w(ops, fin, 1, depth, w, added)
            added[ops[nextop].well][ADDRESS[ops[nextop].reag[0]]] = 1
        else:
            nextop = iddfs_oneiter(ops, fin, 1, depth)

        fin.append(ops[nextop])
        ops.pop(nextop)


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


# single iteration of iddfs (with w)
def iddfs_oneiter_with_w(ops, fin, curdepth, depth, w, added):
    # determine the potential cost of each possible operation
    potcost = []
    for i in range(0, len(ops)):
        potcost.append(cost_func_with_w(fin, ops[i], w, added))
        # next iteration
        if (curdepth < depth and len(ops) != 1):
            pseudoadded = added.copy()
            pseudoadded[ops[i].well][ADDRESS[ops[i].reag[0]]] = 1
            potcost[i] += iddfs_oneiter_with_w(ops[0:i] + ops[(i + 1):len(ops)], fin + [ops[i]], curdepth + 1, depth, w,
                                               pseudoadded)
            # added[ops[i].well][ADDRESS[ops[i].reag[0]]] = 0

    # act according to the determined costs
    if (curdepth == 1):
        answer = potcost.index(min(potcost))
    else:
        answer = min(potcost)

    return answer


# -------------------------------AUXILIARY FUNCTIONS-------------------------------
# get a list of all operations from w
def getops(w, ops):
    for well in range(0, len(w)):
        for reagent in range(0, len(w[well])):
            ops.append(Oper(w[well][reagent], well))


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
