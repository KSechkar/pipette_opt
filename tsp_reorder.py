# PRE-TSP METHOD REORDERINGS
# By Kirill Sechkar
# v0.0.1, 2.7.20

import numpy as np


# ------------------CLASS DEFINITIONS---------------
# needed for the sametogether reordering
class Sametogether:
    def __init__(self, reagtype):
        self.reagtype = reagtype
        self.subs = []
        self.outgoing = 0


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

    for i in range(0, 4):
        print(together[i].outgoing)
