# AUXIL.PY
# Auxiliary functions (e.g. route cost calculator, operation printer, json file reader) used by both TSP- and State Space-based methods
# v0.2.0, 20.8.2020

import numpy as np
from copy import deepcopy


# ------------------------------CLASS DEFINITIONS----------------------------------
# Each part matched with the wells it is added to
class Ss:
    # initialisation
    def __init__(self, part, wellno):
        self.part = part
        self.wells = [wellno]

    # record new well in the subset
    def nuwell(self, wellno):
        self.wells.append(wellno)

    # for printing the subset's part type and wells
    def __str__(self):
        strRep = 'p'+str(self.part[0])+','+str(self.part[1]) + '|'
        for i in range(0, len(self.wells)):
            strRep = strRep + ' ' + str(self.wells[i])
        return strRep


# Final output format is an array of Operations - the same for ALL methods
class Oper:
    # initialisation
    def __init__(self, part, well):
        self.part = part
        self.well = well
        self.changed = False # is True if this operation involves a tip change

    # for printing the part type and destination well
    def __str__(self):
        strRep = 'p'+str(self.part[0])+','+str(self.part[1]) + ' -> w' + str(self.well)
        if(self.changed):
            strRep+=' | change tip'
        return strRep


# ----------------------------------RESULTS PRINTING-------------------------------
def dispoper(fin):
    for i in range(0, len(fin)):
        print(fin[i])


# ------------------------------CAPACITY CALCULATOR--------------------------------
# determine how many doses of each part a pipette can hold (i.e. find capacities)
def capacities(reqvols, pipcap, airgap):
    cap={} # create the capacities dictionary

    # for every part and its required volume, enter the capacity into the dictionary
    for reqvol in reqvols.keys():
        cap[reqvol]=capac(pipcap,reqvols[reqvol],airgap)

    # return the capacities dictionary
    return cap


# capacity calculator for a single part
def capac(pipcap, dose, airgap):
    doseandgap = dose + airgap # calculate how much a dose of part AND the air gap are together
    cap = int(pipcap / doseandgap)  # get how many collections followed by gap will fit
    if (dose <= (pipcap - cap * doseandgap)):  # see if a final collection without an air gap would also fit
        cap += 1
    # return capacity
    return cap


# -------------------------------COST FUNCTIONS------------------------------------
# calculate cost based on tip change indicators of the operations in the sequence
def route_cost(fin):
    cost=0 # initialise the cost variable

    # sum the tip changing indicators (which are 1 if tip is changed, 0 otherwise)
    for op in fin:
        cost+=op.changed

    # return the cost
    return cost


# get route cost independently from the tip change indicators in fin (the well/parts array w is required!)
# use for testing
def validate_cost(fin,w,caps):
    # PART 1: initial preparations

    # make a copy, reset its tip change indicators
    cfin=deepcopy(fin)
    for i in range(0, len(cfin)):
        cfin[i].changed = False

    # array listing operations that have different tip changing status than originally given
    # Note: even if it is non-empty, this does not necessarily mean the original caluclation is wrong
    different=[]

    # create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])),dtype=bool) # added[i][j]==True if part at address j in well i has been added


    # PART 2: get the cost

    # PART 2.1: the first operation in fin
    cost=1 # beginning the distribuiton => new tip taken
    cfin[0].changed = True # indicate the tip's been changed
    added[cfin[0].well][cfin[0].part[0]] = 1 # indicate the part's been added

    # PART 2.2: all other operations
    for i in range(1, len(cfin)):
        one_cost = cost_func_with_w(cfin[0:i], cfin[i], w, added, caps) # get operation cost
        cost += one_cost # add operation cost

        added[cfin[i].well][cfin[i].part[0]] = 1 # indicate the part's been added
        # if the tip's been changed (operation cost 1), indicate that
        if(one_cost==1):
            cfin[i].changed = True

        # if the independently found tip-changing indicator is different, record!
        if (fin[i].changed != cfin[i].changed):
            different.append((cfin[i],str(i))) # record a tuple giving the operation and its place in the sequence

    # PART 3: return the cost and the array of differing operations
    return cost, different


# cost of a single operation op; uses the well array w and the array added giving status of the wells
def cost_func_with_w(fin, op, w, added, caps):
    # if no operations have been performed at all, new tip obiously needed => just return 1
    if(len(fin)==0):
        return 1

    # find the destination well of the last operation performed
    lastwell = fin[-1].well

    # if the next operation involves a different reagent from the last one, new tip needed => cost=1
    if(fin[-1].part!=op.part):
        cost=1

    # otherwise, need a detailed investigation
    else:
        cost = 0

        # check if the last well had any parts the next well does not
        for i in range(0,len(w[0])):
            if(w[lastwell][i]!=fin[-1].part):
                if not (added[lastwell][i]==0 or (w[lastwell][i]==w[op.well][i] and added[op.well][i]==1)):
                    cost=1

        # take into account pipette capacity, IF working with a capacitated problem
        if(caps!=None):
            # only need to do that if cost is ostensibly 0 and the number of operations is less than the capacity
            if((cost==0) and (len(fin)>=caps[op.part])):
                # see how many doses of current part have been delivered
                backforcap = 0
                while (backforcap<len(fin)):
                    backforcap += 1
                    if(fin[-backforcap].changed):
                        break

                # if there is now more room for one more dose for delivery, will have to get a new tip
                if(backforcap==caps[op.part]):
                    cost=1

    return cost


# -----------------------------INPUT CONVERSION------------------------------------
"""
The input can be stored as:
1) a 2d-list w, 
2) list of all operations to do ops,
3) subsets of operations grouped by part (used for the LP method and reorderings)

Conversion between these data types is often necessary
"""

# convert w into subsets
def w_to_subsets(w,subsets):
    for i in range(0, len(w)):
        for j in range(0, len(w[0])):
            # check if there already is a subset for the current reagent
            match = False
            for k in range(0, len(subsets)):
                if (subsets[k].part == w[i][j]):
                    match = True
                    break

            # if yes, just add the current well into the relevant subset
            if (match):
                subsets[k].nuwell(i)
            # if no, create a new subset
            else:
                subsets.append(Ss(w[i][j], i))


# convert subsets into operations list
def subsets_to_ops(subsets,ops):
    for i in range(0, len(subsets)):
        for j in range(0, len(subsets[i].wells)):
            ops.append(Oper(subsets[i].part, subsets[i].wells[j]))


# convert the operations list into subsets
def ops_to_subsets(ops, subsets):
    for op in ops:
        # check if there's a subset for the part already
        match = False
        for subset in subsets:
            if (op.part == subset.part):
                match = True
                break

        # if yes, just add the current well into the relevant subset
        if (match):
            subset.nuwell(op.well)
        # if no, create a new subset
        else:
            subsets.append(Ss(op.part,op.well))
