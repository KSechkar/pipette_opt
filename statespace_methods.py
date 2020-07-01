# STATE SPACE-BASED METHODS (SHORTEST PATH/TREE SEARCH) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.1, 30.6.20

import numpy as np
import time

# import functions from own files
from input_generator import wgenerator
from auxil import route_cost, cost_func, dispoper


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
    ops = []  # an Oper list of operations to be performed
    fin = []  # an Oper list of operations in the order which they should be performed

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to dfine the number of wells
    change 4 last arguments to define the size of p, r, c and t reagent sets"""
    w=wgenerator(96,6,6,3,4)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # use iddfs to solve the problem
    iddfs(w, fin, 2)

    dispoper(fin)
    print('The total number of pipette tips used is ' + str(route_cost(fin)))

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')


# -------------------------------SOLVERS-------------------------------
# iddfs function
def iddfs(w, fin, depth):
    ops = []
    getops(w, ops)

    # if we want to randomise the operation order
    np.random.shuffle(ops)

    all_operations = len(w) * len(w[0])

    fin.append(ops[0])
    ops.pop(0)

    while (len(fin) < all_operations):
        nextop=iddfs_oneiter(ops, fin, 1, depth)
        fin.append(ops[nextop])
        ops.pop(nextop)


# single iteration of iddfs
def iddfs_oneiter(ops, fin, curdepth, depth):
    # determine the potential cost of each possible operation
    potcost = []
    for i in range(0, len(ops)):
        potcost.append(cost_func(fin, ops[i]))
        # next iteration
        if (curdepth < depth and len(ops)!=1):
            pseudofin = fin.copy()
            pseudofin.append(ops[i])
            pseudoops = ops.copy()
            pseudoops.pop(i)
            potcost[i] += iddfs_oneiter(pseudoops, pseudofin, curdepth + 1, depth)

    # act according to the determined costs
    if(curdepth==1):
        answer=potcost.index(min(potcost))
    else:
        answer=min(potcost)
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
