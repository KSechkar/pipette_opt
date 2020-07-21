# PRE-TSP METHOD REORDERINGS
# By Kirill Sechkar
# v0.1.0.lp, 15.7.20

import numpy as np
from tsp_lp_solver import tsp_lp_gurobi

from tspy import TSP
from tspy.solvers.utils import get_cost


# ------------------CLASS DEFINITIONS---------------
# needed for the sametogether reordering
class Sametogether:
    def __init__(self, reagtype):
        self.reagtype = reagtype
        self.subs = []
        self.outgoing = 0


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


# ------------------MAIN (only for function testing)---------------
def main():
    D = np.zeros((3, 3))
    subset = Ss('p1', 1)
    subset.wells.append(2)
    print(Dupdate(D, subset))
    print(D)


# -------------SIMPLE REORDERINGS----------------
# totalwells is the total number of wells in the input
def leastout(subsets, totalwells):
    # determine the number of outgoing edges for each subset
    for i in range(0, len(subsets)):
        subsets[i].outgoing = len(subsets[i].wells) * (totalwells - len(subsets[i].wells))

    # sort the list
    subsets.sort(key=lambda subsets: subsets.outgoing)


def sametogether(subsets, totalwells):
    # intialise the 4 lists of same-type reagent subsets
    together = [Sametogether('p'), Sametogether('r'), Sametogether('c'), Sametogether('t')]

    # distribute the subsets among the lists, calculate the total number of outgoing edges for each list
    for i in range(0, len(subsets)):
        # find into which of 4 list we put it
        for position in range(0, 4):
            if (subsets[i].reag[0] == together[position].reagtype):
                break

        together[position].subs.append(subsets[i])  # record the subset in the proper array
        together[position].outgoing += len(subsets[i].wells) * (totalwells - len(
            subsets[i].wells))  # update number of outgoing edges (for further OPTIONAL sorting)

    # sort the 4 lists by the number of outgoing edges (OPTIONAL)
    together.sort(key=lambda together: together.outgoing)

    # record the rearranged subsets
    inlist = 0  # counter within one of the 4 lists
    whichlist = 0  # which of the 4 lists is current
    for i in range(0, len(subsets)):
        subsets[i] = together[whichlist].subs[inlist]
        inlist += 1
        if (inlist == len(together[whichlist].subs)):
            whichlist += 1
            inlist = 0


# -------------STATE-SPACE REORDERINGS----------------
# iddfs
def reorder_iddfs(origsubs, subsets, D, depth):
    # if we want to randomise the order first
    # np.random.shuffle(origsubs)

    all_operations = len(origsubs)

    subsets.append(origsubs[0])
    origsubs.pop(0)
    Dupdate(D, subsets[-1])

    while (len(subsets) < all_operations):
        #print(len(subsets))

        nextop = reorder_iddfs_oneiter(origsubs, subsets, D.copy(), 1, depth)
        subsets.append(origsubs[nextop])
        origsubs.pop(nextop)
        Dupdate(D, subsets[-1])


def reorder_iddfs_oneiter(origsubs, subsets, D, curdepth, depth):
    # determine the potential cost of each possible operation
    potcost = []
    for i in range(0, len(origsubs)):
        potcost.append(solveforcost(origsubs[i], D))
        # next iteration
        werenew = Dupdate(D, origsubs[i])
        if (curdepth < depth and len(origsubs) != 1):
            # change the inputs for next iteration
            subsets.append(origsubs[i])
            origsubs.pop(i)

            # call next iteration
            potcost[i] += reorder_iddfs_oneiter(origsubs, subsets, D, curdepth + 1, depth)

            # change the inputs back
            origsubs.insert(i, subsets[-1])
            subsets.pop()
        Drollback(D, werenew)

    # act according to the determined costs
    if (curdepth == 1):
        answer = potcost.index(min(potcost))
    else:
        answer = min(potcost)

    return answer


# greedy algorithm
def reorder_greedy(origsubs, subsets, D, heur):
    # if we want to randomise the order first
    # np.random.shuffle(origsubs)

    all_operations = len(origsubs)

    subsets.append(origsubs[0])
    origsubs.pop(0)
    Dupdate(D, subsets[0])

    while (len(subsets) < all_operations):
        #print(len(subsets))

        nextop = reorder_greedy_onestep(origsubs, subsets, D, heur)
        subsets.append(origsubs[nextop])
        origsubs.pop(nextop)
        Dupdate(D, subsets[-1])


def reorder_greedy_onestep(origsubs, subsets, D, heur):
    # determine the potential cost of each possible operation
    potcost = []
    for i in range(0, len(origsubs)):
        # cost function component
        potcost.append(solveforcost(origsubs[i], D))

        # heuristic component
        subsets.append(origsubs[i])
        origsubs.pop(i)
        #print(h_tree(subsets, origsubs, D.copy(), heur))  # TEST ONLY
        potcost[-1] += h_tree(subsets, origsubs, D.copy(), heur)
        origsubs.insert(i, subsets[-1])
        subsets.pop()

    return potcost.index(min(potcost))


# a star
def reorder_a_star(origsubs, subsets, D, heur):
    alloperations = len(origsubs)

    # if we want to randomise the operation order
    # np.random.shuffle(ops)

    # starting position - state with no operations performed
    states = [[]]  # list of states in state-space under consideration
    unstates = [origsubs]  # list of reagents NOT added for a given state
    g = [0]  # distance from origin, mirrors states
    h = [0]  # heuristic function values
    f = [0]  # f(states[i])=g(states[i])+h(states[i]), mirrors states
    l = [0]
    d = [D.copy()]

    while (len(states) != 0):
        consider = f.index(min(f))

        # TEST ONLY
        #if (len(states[consider]) != 0):
            #print(str(len(states)) + ' ' + str(h[consider]) + ' ' + str(max(l)))

        #print(consider)
        if (len(states[consider]) == alloperations):
            for s in states[consider]:
                subsets.append(s)
            break

        # add neighbours to considered states
        for i in range(0, len(unstates[consider])):
            states.append(states[consider] + [unstates[consider][i]])
            d.append(D.copy())
            unstates.append(unstates[consider][0:i] + unstates[consider][i + 1:len(unstates[consider])])
            g.append(g[consider] + solveforcost(unstates[consider][i], d[-1]))
            h.append(h_tree(states[-1], unstates[-1], d[-1], heur))
            f.append(g[-1] + h[-1])
            Dupdate(d[-1], states[-1][-1])
            l.append(len(states[-1]))  # TEST ONLY

        # remove the state in question
        states.pop(consider)
        unstates.pop(consider)
        g.pop(consider)
        f.pop(consider)
        d.pop(consider)
        l.pop(consider)  # TEST ONLY


def h_tree(complete, todo, D, heur):
    wt = 1  # weighing factor on the whole heuristic
    if (heur == 'countall'):  # h is the sum of all nonzero edges multiplied by the number of iterations left
        k = 1  # weighting factor of tally vs the amount left
        Dupdate(D, complete[-1])
        tally = 0
        for i in range(0, len(D)):
            for j in range(0, len(D)):
                if (i != j):
                    tally += D[i][j]
        answer = wt * ((tally ** k) * len(todo))# ** (1 / (k + 1))
        return answer
    elif (heur == 'leastout'):  # works analogically to the leastout sorting
        k = 2
        l = len(Dupdate(D, complete[-1]))
        answer = wt * ((l ** k) * len(todo)) ** (1 / (k + 1))
        return answer

    return 0


# -----------STATE-SPACE auxiliaries--------------
# update D according to the subset
def Dupdate(D, subset):
    werenew = []
    sublen = len(subset.wells)
    for i_well in range(0, sublen):
        current_well = 0
        for j_D in range(0, len(D)):
            if (j_D == subset.wells[current_well]):
                if (current_well < sublen - 1):
                    current_well += 1
            else:
                if (D[subset.wells[i_well]][j_D] != 1):
                    D[subset.wells[i_well]][
                        j_D] = 1  # make the edge going from the subset into the rest of D equal to one (updating D)
                    werenew.append([subset.wells[i_well], j_D])
    return werenew


# roll D back to before the Dupdate call
def Drollback(D, werenew):
    for i in range(0, len(werenew)):
        D[werenew[i][0]][werenew[i][1]] -= 1


# obtain a solution for given subset to determine its cost for sure
def solveforcost(subset, D):
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
                subD[i_well + 1][current_well + 1] = D[subset.wells[i_well]][j_D]  # select the edges within the subset into the submatrix
                if (current_well < sublen - 1):
                    current_well += 1

    # PART 3: solve TSP for the subset
    tour=tsp_lp_gurobi(subD)

    # PART 4: record the operations into the final output, 'unwrapping' the cycle arround the added zero node to create a path
    # find the position of the zero node in the tour
    for i in range(1,len(tour)):
        fin.append(Oper(subset.reag, subset.wells[tour[i] - 1]))

    # PART 5: return the adjusted number of pipette tip changes
    tour.append(0) # accounts for notation difference between tspy and gurobi
    tsp = TSP()
    tsp.read_mat(subD)
    return tips + get_cost(tour, tsp)  # include the tour cost in the number of tip changes


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
