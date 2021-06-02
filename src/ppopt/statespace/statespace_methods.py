# STATE SPACE-BASED METHODS (SHORTEST PATH/TREE SEARCH) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.1.0, 22.7.20

from ppopt.auxil import *


# ---------------------------------NNs SOLVER------------------------------------
# Nearest Neighbour tree search; the depth argument determines search depth
def nns(w, fin, depth, reord,caps):
    # PART 1: intial preparations

    # PART 1.1: get an Oper list of operations to be performed
    ops = []
    w_to_ops(w, ops, reord)
    all_operations = len(ops) # get the total number of subsets

    # PART 1.3: make w and capacity global for simplicity
    global globw, globcaps
    globw = w
    globcaps = caps

    # PART 1.4: create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])), dtype=bool)  # added[i][j]==True if part at address j in well i has been added


    # PART 2: get the sequence of operations

    # PART 2.1: first operation
    fin.append(ops[0]) # record operation
    added[fin[0].well][fin[0].part[0]] = 1  # indicate the part's been added
    fin[0].changed = True # beginning the distribuiton => new tip taken => change the indicator
    ops.pop(0) # remove from the list of unperformed operations

    # PART 2.2: all other operations
    while (len(fin) < all_operations):
        #print(str(len(fin))+' of '+str(all_operations)+' operations')  # uncomment if need to track the progress
        # get next operation
        nextop = nns_oneiter_with_w(ops, fin, 1, depth,added)
        nextcost = cost_func_with_w(fin, ops[nextop], w, added, caps) # get next operation's cost

        fin.append(ops[nextop]) # record next operation
        added[ops[nextop].well][ops[nextop].part[0]] = 1 #indicate the part's been added
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
            added[ops[i].well][ops[i].part[0]] = 1
            fin.append(ops[i])
            fin[-1].changed=True
            ops.pop(i)

            # call the single-iteration function again
            potcost[i] += nns_oneiter_with_w(ops, fin, curdepth + 1, depth, added)

            # change the inputs back
            ops.insert(i, fin[-1])
            fin[-1].changed = False
            fin.pop()
            added[ops[i].well][ops[i].part[0]] = 0

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

    # PART 1.1: get an Oper list of operations to be performed
    ops = []
    w_to_ops(w, ops, reord)
    all_operations = len(ops)  # get the total number of subsets

    # PART 1.3: make w and capacity global for simplicity
    global globw, globcaps
    globw = w
    globcaps = caps

    # PART 1.4: create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])), dtype=bool)  # added[i][j]==True if part at address j in well i has been added


    # PART 2: get the sequence of operations

    # PART 2.1: first operation
    fin.append(ops[0])  # record operation
    added[fin[0].well][fin[0].part[0]] = 1  # indicate the part's been added
    fin[0].changed = True  # beginning the distribuiton => new tip taken => change the indicator
    ops.pop(0)  # remove from the list of unperformed operations

    # PART 2.2: all other operations
    while (len(fin) < all_operations):
        print(str(len(fin))+' of '+str(all_operations)+' operations')  # uncomment if need to track the progress

        # get next operation
        nextop = greedy_tree_onestep(ops, fin, w, added, heur)
        nextcost = cost_func_with_w(fin, ops[nextop], w, added, caps)  # get next operation's cost

        fin.append(ops[nextop])  # record next operation
        added[ops[nextop].well][ops[nextop].part[0]] = 1  # indicate the part's been added
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
            # for the current part (subset), get the number of tips needed to deliver
            est_cost += np.ceil(len(subset.wells)/globcaps[subset.part])

            # if the current part (subset) is the same as last part added, assume the tip is kept from then
            if(subset.part == state[-1].part):
                est_cost -= 1
        return est_cost


# -------------------------------MAIN (TESTING ONLY!)-----------------------------
def main():
    fin = []  # final array where the operations are to be recorded

    """
    INPUT:
    a) Use a manually defined well array, or
    b) Generate 
            change 1st argument of wgenerator to define the number of wells/constructs
            change 4 last arguments of wgenerator to define the size of p, r, c and t part sets (this is Start-Stop assembly)
    Comment out the respective line to deselect
    """
    w = [[(0, 1), (1, 2), (2, 4), (3, 1)],
         [(0, 2), (1, 2), (2, 1), (3, 1)],
         [(0, 1), (1, 2), (2, 2), (3, 2)],
         [(0, 2), (1, 3), (2, 1), (3, 1)]]

    w = wgenerator(96, 6, 6, 3, 4)

    # generate required volumes (for testing). Values taken from a real instance of Start-Stop assembly
    ss = []
    w_to_subsets(w, ss)
    reqvols = {}
    for s in ss:
        if (s.part[0] == 0):
            reqvols[s.part] = 1.09
        elif (s.part[0] == 1):
            reqvols[s.part] = 0.33
        elif (s.part[0] == 2):
            reqvols[s.part] = 0.36
        else:
            reqvols[s.part] = 0.75

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # get capacitites
    caps = capacities(reqvols, 10, 1.0)

    # DECOMMENT THE LINE THAT CALLS THE ALGORITHM YOU WANT TO TEST, COMMENT ALL OTHERS
    # CHANGE THE 'REORD' ARGUMENT TO APPLY AN HEURSTIC REORDERING BEFORE SOLVING...
    # ...reord=None FOR NO REORDERING; reord='sameotogether' FOR SAMETOGETHER REORDERING.

    # use nearest-neighbour tree search (NNS) to solve the problem
    # nns(w, fin, 1, reord='sametogether', caps=caps)

    # use depth 2 NNS to solve the problem
    # nns(w, fin, 2, reord='sametogether', caps=caps)

    # use greedy algorithm on a tree to solve the problem
    greedy_tree(w, fin, 'optimistic+cap', reord='sametogether', caps=caps)

    # display the solution
    dispoper(fin)

    # PERFORMANCE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')
    print('The total number of pipette tips used is (from resultant list) ' + str(route_cost(fin)))
    print('The total number of pipette tips used is (independent caluclation) ' + str(independent_cost(fin, w, caps)[0]))


if __name__ == "__main__":
    # import necessary modules for independent test running
    import time
    from src.ppopt import wgenerator
    main()
