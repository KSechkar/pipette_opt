# AUXIL.PY
# Auxiliary functions (e.g. route cost calculator, operation printer) used by both TSP- and State Space-based methods
# v0.0.1, 01.07.2020

import numpy as np

# -------------RESULTS PRINTING-------------
def dispoper(fin):
    for i in range(0, len(fin)):
        print(fin[i])


# ---------------------COSTS----------------
# return route cost
def route_cost(fin):
    cost = 1
    for f in range(1, len(fin)):
        cost += cost_func(fin[0:f], fin[f])
    return cost


# return cost of an operation in the light of all previous operations (fin)
def cost_func(fin, op):
    # find the index of last for easier further referencing
    lastindex = len(fin) - 1
    # find the number of the last well where a reagent was added
    lastwell = fin[lastindex].well
    # reagents present in this well affect the cost and thus need to be determined
    lastwell_addedreags = []
    for i in range(0, lastindex):
        if (fin[i].well == lastwell):
            lastwell_addedreags.append(fin[i].reag)

    # determine the cost of each possible operation
    cost = 1  # by default, it is 1 for any operation
    if (op.reag == fin[lastindex].reag):
        # this part of the program checks the conditions outlined in the 'Problem Representation' Document
        for i in range(0, lastindex):
            for j in range(0, len(lastwell_addedreags)):
                if (fin[i].well == op.well) and (fin[i].reag == lastwell_addedreags[j]):
                    lastwell_addedreags.pop(j)
                    break
            if (len(lastwell_addedreags) == 0):
                break
        if (len(lastwell_addedreags) == 0):
            cost = 0
    return cost

#INCOMPLETE
# a faster way to determine operation cost by also using the input list and knowing which position in w stands for what
def cost_func_with_w(w,fin,op):
    reag_address={'p': 0, 'r': 1, 'c': 2, 't': 3}
    # find the index of last for easier further referencing
    lastindex = len(fin) - 1
    # find the number of the last well where a reagent was added
    lastwell = fin[lastindex].well

    #generate a table of whether each reagent in each well was added
    added=np.zeros(len(w),len(w[0]))
    for f in fin:
        added[f.well][reag_address[f.reag[0]]]=1

