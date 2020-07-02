# AUXIL.PY
# Auxiliary functions (e.g. route cost calculator, operation printer) used by both TSP- and State Space-based methods
# v0.0.2, 02.07.2020

import numpy as np

#match reagents with their addresses in w
ADDRESS={'p': 0, 'r': 1, 'c': 2, 't': 3}

# -------------RESULTS PRINTING-------------
def dispoper(fin):
    for i in range(0, len(fin)):
        print(fin[i])


# ---------------------COSTS----------------
# return route cost
def route_cost(fin):
    cost = 1
    for f in range(1, len(fin)):
        cost+=cost_func(fin[0:f], fin[f])
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


# alterative way to determine operation cost by also knowing the the input well array
def route_cost_with_w(fin,w):
    added = np.zeros((len(w), len(w[0]))) #tells which reagents were added to which well

    #for the first operation in fin
    cost=1
    added[fin[0].well][ADDRESS[fin[0].reag[0]]]=1
    for i in range(1, len(fin)):
        cost+=cost_func_with_w(fin[0:i],fin[i],w,added)
        added[fin[i].well][ADDRESS[fin[i].reag[0]]] = 1
    return cost

def cost_func_with_w(fin,op,w,added):
    # find the index of last for easier further referencing
    lastindex = len(fin) - 1
    # find the number of the last well where a reagent was added
    lastwell = fin[lastindex].well

    if(fin[lastindex].reag!=op.reag):
        cost=1
    else:
        #check  on all 3 remaining
        cost=0
        for i in range(0,len(w[0])):
            if(w[lastwell][i]!=fin[lastindex].reag):
                if not (added[lastwell][i]==0 or (w[lastwell][i]==w[op.well][i] and added[op.well][i]==1)):
                    cost=1

    return cost



