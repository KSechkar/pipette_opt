# STATE SPACE-BASED METHODS (SHORTEST PATH/TREE SEARCH) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.1.0, 22.7.20

import time

# import functions from own files
from input_generator import wgenerator
from auxil import *
from tsp_reorder import sametogether, leastout


# -------------------------------CLASS DEFINITIONS-------------------------------
# Final output format is an array of Operations - the same for ALL methods
class Oper:
    # initialisation
    def __init__(self, part, well):
        self.part = part
        self.well = well
        self.changed = False # is True if this operation involves a tip change

    # for printing the part type and destination well
    def __str__(self):
        strRep = self.part + ' -> w' + str(self.well)
        return strRep


# ---------------------------------NNs SOLVER------------------------------------
# Nearest Neighbour tree search; the depth argument determines search depth
def nns(w, fin, depth, reord,caps):
    # PART 1: intial preparations

    # PART 1.1: get part type addresses from w
    address = addrfromw(w)

    # PART 1.2: get an Oper list of operations to be performed
    ops = []
    getops(w, ops, reord)
    all_operations = len(ops) # get the total number of subsets

    # PART 1.3: make w, address and capacity global for simplicity
    global globw, globcaps, globaddress
    globw = w
    globaddress = address
    globcaps = caps

    # PART 1.4: create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])), dtype=bool)  # added[i][j]==True if part at address j in well i has been added


    # PART 2: get the sequence of operations

    # PART 2.1: first operation
    fin.append(ops[0]) # record operation
    added[fin[0].well][address[fin[0].part[0]]] = 1  # indicate the part's been added
    fin[0].changed = True # beginning the distribuiton => new tip taken => change the indicator
    ops.pop(0) # remove from the list of unperformed operations

    # PART 2.2: all other operations
    while (len(fin) < all_operations):
        # get next operation
        nextop = nns_oneiter_with_w(ops, fin, 1, depth,added)
        nextcost = cost_func_with_w(fin, ops[nextop], w, added, caps) # get next operation's cost

        fin.append(ops[nextop]) # record next operation
        added[ops[nextop].well][address[ops[nextop].part[0]]] = 1 #indicate the part's been added
        # if the tip's been changed (operation cost 1), indicate that
        if (nextcost == 1):
            fin[-1].changed = True
        ops.pop(nextop) # remove the operation from the list of unperformed operations


# single iteration of NNs
def nns_oneiter_with_w(ops, fin, curdepth, depth, added):
    # PART 1: determine the potential cost of each possible operation
    potcost = [] # initialise list of potential costs
    for i in range(0, len(ops)):
        # get cost of choosing this operation
        potcost.append(cost_func_with_w(fin, ops[i], globw, added, globcaps))
        # go deeper (if needed and possible)
        if (curdepth < depth and len(ops) != 1):
            # change the inputs as if this operation was chosen
            added[ops[i].well][globaddress[ops[i].part[0]]] = 1
            fin.append(ops[i])
            fin[-1].changed=True
            ops.pop(i)

            # call the single-iteration function again
            potcost[i] += nns_oneiter_with_w(ops, fin, curdepth + 1, depth, added)

            # change the inputs back
            ops.insert(i, fin[-1])
            fin[-1].changed = False
            fin.pop()
            added[ops[i].well][globaddress[ops[i].part[0]]] = 0

    # PART 2: act according to the determined costs
    # if the current depth is 1, return the entry with the least potential cost
    if (curdepth == 1):
        answer = potcost.index(min(potcost))
    # if current depth is greater, return the minimum potential cost to add to current subset's potential cost
    else:
        answer = min(potcost)

    return answer


# ----------------------------------GREEDY SOLVER--------------------------------
# Greedy search; the heur argument determines which heuristic is used
def greedy_tree(w, fin, heur, reord,caps):
    # PART 1: intial preparations

    # PART 1.1: get part type addresses from w
    address = addrfromw(w)

    # PART 1.2: get an Oper list of operations to be performed
    ops = []
    getops(w, ops, reord)
    all_operations = len(ops)  # get the total number of subsets

    # PART 1.3: make w, address and capacity global for simplicity
    global globw, globcaps, globaddress
    globw = w
    globaddress = address
    globcaps = caps

    # PART 1.4: create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])), dtype=bool)  # added[i][j]==True if part at address j in well i has been added


    # PART 2: get the sequence of operations

    # PART 2.1: first operation
    fin.append(ops[0])  # record operation
    added[fin[0].well][address[fin[0].part[0]]] = 1  # indicate the part's been added
    fin[0].changed = True  # beginning the distribuiton => new tip taken => change the indicator
    ops.pop(0)  # remove from the list of unperformed operations

    # PART 2.2: all other operations
    while (len(fin) < all_operations):
        # get next operation
        nextop = greedy_tree_onestep(ops, fin, w, added, heur)
        nextcost = cost_func_with_w(fin, ops[nextop], w, added, caps)  # get next operation's cost

        fin.append(ops[nextop])  # record next operation
        added[ops[nextop].well][address[ops[nextop].part[0]]] = 1  # indicate the part's been added
        # if the tip's been changed (operation cost 1), indicate that
        if (nextcost == 1):
            fin[-1].changed = True
        ops.pop(nextop)  # remove the operation from the list of unperformed operations


# single step of greedy search
def greedy_tree_onestep(ops, fin, w, added, heur):
    # PART 1: determine the potential cost+heuristic of each possible operation
    potcost = []
    for i in range(0, len(ops)):
        # PART 1.1: cost function component
        potcost.append(cost_func_with_w(fin, ops[i], w, added, globcaps))

        # PART 1.2: heuristic component
        # pretend the current subset is included in the route
        fin.append(ops[i])
        ops.pop(i)

        # get the heuristic value
        potcost[-1] += h_tree(fin, ops, heur)

        # undo including the current subset into the route
        ops.insert(i, fin[-1])
        fin.pop()

    # PART 2: return the index of subset with the minimal cost+heuristic
    return potcost.index(min(potcost))


# heuristic function
def h_tree(state, unstate, heur):
    # optimistic: h is the number of future tip changes assuming tips only change between different parts
    if (heur == 'optimistic'):
        already = []
        for unop in unstate:
            ispresent = False
            for alread in already:
                if (unop.part == alread):
                    ispresent = True
                    break
            if not ispresent:
                already.append(unop.part)
        return len(already)

    # optimistic+cap: h is the number of future  tip changes,
    # assuming tips only change between different parts and due to pipette capacity
    elif (heur == 'optimistic+cap'):
        subsets=[]
        ops_to_subsets(unstate,subsets)
        est_cost=0
        for subset in subsets:
            est_cost+=1
            est_cost+=round(len(subset.wells)/globcaps[subset.part])
        return est_cost


# ------------------------------AUXILIARY FUNCTIONS-------------------------------
# get a list of all operations from w
def getops(w, ops, reord):
    # no reordering => just convert w into ops
    if (reord == None):
        for well in range(0, len(w)):
            for part in range(0, len(w[well])):
                ops.append(Oper(w[well][part], well))
        return

    # random reordering => convert into ops, randomly shuffle ops
    elif (reord == 'random'):
        for well in range(0, len(w)):
            for part in range(0, len(w[well])):
                ops.append(Oper(w[well][part], well))
        np.random.shuffle(ops)
        return

    # subset-based reorderings => convert into subsets, apply a corresponding reordering, convert into ops
    if (reord == 'leastout' or reord == 'sametogether' or reord == 'justsubsets'):
        subsets = []
        w_to_subsets(w, subsets)
        if (reord == 'sametogether'):
            sametogether(subsets, w)
        else:
            leastout(subsets,w)
        subsets_to_ops(subsets, ops)


# -------------------------------MAIN (TESTING ONLY!)-----------------------------
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

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # generate required volumes (for testing)
    ss = []
    w_to_subsets(w, ss)
    reqvols = {}
    for s in ss:
        if (s.part[0] == 'p'):
            reqvols[s.part] = 1.09
        elif (s.part[0] == 'r'):
            reqvols[s.part] = 0.33
        elif (s.part[0] == 'c'):
            reqvols[s.part] = 0.36
        else:
            reqvols[s.part] = 0.75

    # get capacitites
    caps = capacities(reqvols, 10, 1.0)

    # use nearest-neighbour tree search (NNS) to solve the problem
    # nns(w, fin, 1, reord='sametogether', caps=caps)

    # use depth 2 NNS to solve the problem
    nns(w, fin, 2, reord='sametogether', caps=caps)

    # use greedy algorithm on a tree to solve the problem
    # greedy_tree(w, fin, 'optimistic+cap', reord='sametogether', caps=caps)

    dispoper(fin)

    # PERFORMANCE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')

    print('The total number of pipette tips used is ' + str(route_cost_with_w(fin, w, caps)))


if __name__ == "__main__":
    main()
