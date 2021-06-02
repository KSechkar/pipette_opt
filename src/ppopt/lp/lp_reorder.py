# HEURISTIC REORDERINGS OF THE DNA PART LIST
# these are the reorderings only for the LP-based method, which are based on tree search algorithms
# By Kirill Sechkar
# v0.1.0, 30.5.21

import numpy as np
from tspy import TSP
from tspy.solvers.utils import get_cost
from tspy.solvers import TwoOpt_solver

from src.ppopt import lp_cap

# Nearest Neighbour tree search; the depth argument determines search depth
def reorder_nns(origsubs, subsets, D, depth,caps):
    # PART 1: initial preparations
    # set maximum optimisation time for the cost function (less than default for quicker searching)
    global maxtime
    maxtime=0.1

    # get the total number of subsets
    all_subsets = len(origsubs)


    # PART 2: reorder the subsets

    # PART 2.1: the first subset
    subsets.append(origsubs[0])
    origsubs.pop(0)
    Dupdate(D, subsets[-1])

    # PART 2.2: all other subsets
    while (len(subsets) < all_subsets):
        # get next subset
        nextsub = reorder_nns_oneiter(origsubs, subsets, D.copy(), 1, depth, caps)

        subsets.append(origsubs[nextsub]) # record next subset
        origsubs.pop(nextsub) # remove next subset from the list of unvisited subsets
        Dupdate(D, subsets[-1]) # update D


# single iteration of NNs
def reorder_nns_oneiter(origsubs, subsets, D, curdepth, depth, caps):
    # PART 1: determine the potential cost of each possible operation
    potcost = [] # initialise list of potential costs
    for i in range(0, len(origsubs)):
        # get cost of choosing this susbset
        potcost.append(solveforcost(origsubs[i], D,caps[origsubs[i].part]))
        # go deeper (if needed and possible)
        if (curdepth < depth and len(origsubs) != 1):
            # change D and inputs as if the current subset was chosen
            werenew = Dupdate(D, origsubs[i])
            subsets.append(origsubs[i])
            origsubs.pop(i)

            # call the single-iteration function again
            potcost[i] += reorder_nns_oneiter(origsubs, subsets, D, curdepth + 1, depth,caps)

            # change the inputs back
            origsubs.insert(i, subsets[-1])
            subsets.pop()
            Drollback(D, werenew)

    # PART 2: act according to the determined costs
    # if the current depth is 1, return the entry with the least potential cost
    if (curdepth == 1):
        answer = potcost.index(min(potcost))
    # if current depth is greater, return the minimum potential cost to add to current subset's potential cost
    else:
        answer = min(potcost)

    return answer


# Greedy search; the heur argument determines which heuristic is used
def reorder_greedy(origsubs, subsets, D, heur, caps):
    # PART 1: initial preparations
    # set maximum optimisation time for the cost function  (less than default for quicker searching)
    global maxtime
    maxtime = 0.1

    # get the total number of subsets
    all_subsets = len(origsubs)


    # PART 2: reorder the subsets

    # PART 2.1: the first subset
    subsets.append(origsubs[0])
    origsubs.pop(0)
    Dupdate(D, subsets[0])

    # PART 2.2: all other subsets
    while (len(subsets) < all_subsets):
        # get next subset
        nextop = reorder_greedy_onestep(origsubs, subsets, D, heur,caps)

        subsets.append(origsubs[nextop]) # record next subset
        origsubs.pop(nextop) # remove next subset from the list of unvisited subsets
        Dupdate(D, subsets[-1]) # update D


# single step of greedy search
def reorder_greedy_onestep(origsubs, subsets, D, heur,caps):
    # PART 1: determine the potential cost+heuristic of each possible operation
    potcost = [] # intialise the list of potetntial values
    for i in range(0, len(origsubs)):
        # PART 1.1: cost function component
        potcost.append(solveforcost(origsubs[i], D,caps[origsubs[i].part]))

        # PART 1.2: heuristic component
        # pretend the current subset is included in the route (D NOT updated! As some heuristics MIGHT need unupdated D)
        subsets.append(origsubs[i])
        origsubs.pop(i)

        # get the heuristic value
        potcost[-1] += h_tree(subsets, origsubs, D.copy(), heur)

        # undo including the current subset into the route
        origsubs.insert(i, subsets[-1])
        subsets.pop()

    # PART 2: return the index of subset with the minimal cost+heuristic
    return potcost.index(min(potcost))

# heursiric function
# complete = subsets currently in route, remaining = subsets not in route
def h_tree(inroute, remaining, D, heur):
    # countall: h is the sum of all edge costs multiplied by the number of iterations left
    if (heur == 'countall'):
        Dupdate(D, inroute[-1]) # update D for current subset (as it was not updated before)

        # sum all edge costs
        tally = 0
        for i in range(0, len(D)):
            for j in range(0, len(D)):
                if (i != j):
                    tally += D[i][j]

        # calculate the value of h
        answer = tally * len(remaining)
        return answer


# ----------------------------AUXILIARIES FOR TREE SEARCH---------------------------
# update D according to the subset
def Dupdate(D, subset):
    werenew = [] # created the list of altered edges
    sublen = len(subset.wells) # get subgraph length to avoid calling len too often

    """
    For every well belonging to the subset (given by subset.wells[i_well]), we consider all outgoing edges
    (given by D[subset.wells[i_well]][j_D]).
    If the edge arrives to a well not in the subset, its cost in D must be updated to be 1 (if it wasn't 1 already).
    """
    for i_well in range(0, sublen):
        current_well = 0
        for j_D in range(0, len(D)):
            if (j_D == subset.wells[current_well]):
                if (current_well < sublen - 1):
                    current_well += 1
            else:
                if (D[subset.wells[i_well]][j_D] != 1):
                    D[subset.wells[i_well]][j_D] = 1
                    werenew.append([subset.wells[i_well], j_D]) # record which edge's cost was altered

    # return the list of all altered edges
    return werenew


# roll D back to before the Dupdate call
def Drollback(D, werenew):
    for i in range(0, len(werenew)):
        D[werenew[i][0]][werenew[i][1]] -= 1


# obtain a solution for given subset to determine its cost for sure
def solveforcost(subset, D, cap):
    # PART 1: initial preparations
    # get length to avoid calling len too often
    sublen = len(subset.wells)

    # initialise the subset's matrix subD
    subD = np.zeros((sublen + 1, sublen + 1))  # an extra 0 node is needed for tspy compatibility

    # PART 2: select submatrix and update D as if problem for the subset is already solved
    """
    For every well belonging to the subset (given by subset.wells[i_well]), we consider all outgoing edges
    (given by D[subset.wells[i_well]][j_D]).
    If the edge arrives at another subset well (given by subset.wells[current_well]), we select it into the subgraph,
    i.e. we write its cost into subD[i_well+1][current_well].
    """
    for i_well in range(0, sublen):
        current_well = 0
        for j_D in range(0, len(D)):
            if (j_D == subset.wells[current_well]):
                subD[i_well + 1][current_well + 1] = D[subset.wells[i_well]][j_D]
                if (current_well < sublen - 1):
                    current_well += 1

    # PART 3: solve TSP for the subset and record costs
    # 3a): capacitated problem
    if (cap!=None):
        # get the chain coverage
        if (len(subD) == 2):
            chains = [[1]]
        else:
            chains = lp_cap(subD, cap, maxtime)

        # the number of chains is the number of tips used
        cost = len(chains)

    # 3b): non-capacitated problem, uses the tspy package
    else:
        # get the TSP tour using the tspy package
        tsp = TSP()
        tsp.read_mat(subD)
        two_opt = TwoOpt_solver(initial_tour='NN', iter_num=100)
        tour = tsp.get_approx_solution(two_opt)

        # get cost of the tour
        cost = get_cost(tour, tsp)

    # PART 4: return the cost
    return cost


# ------------------------------MAIN (TESTING ONLY)--------------------------------
def main():
    D = np.zeros((3, 3))
    subset = Ss('p1', 1)
    subset.wells.append(2)
    print(Dupdate(D, subset))
    print(D)


if __name__ == "__main__":
    from src.ppopt import Ss
    main()
