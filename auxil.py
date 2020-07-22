# AUXIL.PY
# Auxiliary functions (e.g. route cost calculator, operation printer, json file reader) used by both TSP- and State Space-based methods
# v0.1.1, 22.7.2020

import numpy as np
import json

#match reagents with their addresses in w
ADDRESS={'p': 0, 'r': 1, 'c': 2, 't': 3}

# -----------------------class definitions-----------------------------
#operations
class Oper:
    def __init__(self, reag, well):
        self.reag = reag
        self.well = well

    def __str__(self):  # for printing the subset's reagent type and wells out
        strRep = self.reag + ' -> w' + str(self.well)
        return strRep

#operation subsets
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


# -----------------------RESULTS PRINTING-----------------------------
def dispoper(fin):
    for i in range(0, len(fin)):
        print(fin[i])


# -----------------------RESULTS PRINTING-----------------------------
# need to clairfy if the last collection has an air gap after it
# currently assumed that NO
def capac(pipcap, dose, airgap):
    doseandgap = dose + airgap
    cap = int(pipcap / doseandgap)  # get how many collections followed by gap will fit
    if (dose <= (pipcap - cap * doseandgap)):  # see if a final collection without air gap would also fit
        cap += 1
    return cap


# ---------------------------COST FUNCTIONS------------------------------------
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

    if (lastindex < 0):  # if no previous operations have been performed, we obviously need to put on a tip
        return 1

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
def route_cost_with_w(fin,w,cap):
    added = np.zeros((len(w), len(w[0]))) #tells which reagents were added to which well

    #for the first operation in fin
    cost=1
    added[fin[0].well][ADDRESS[fin[0].reag[0]]]=1
    for i in range(1, len(fin)):
        cost += cost_func_with_w(fin[0:i], fin[i], w, added,cap)
        added[fin[i].well][ADDRESS[fin[i].reag[0]]] = 1
    return cost


def cost_func_with_w(fin,op,w,added,cap):
    # find the index of last for easier further referencing
    lastindex = len(fin) - 1
    if(lastindex<0): #if no previous operations have been performed, we obviously need to put on a tip
        return 1

    # find the num ber of the last well where a reagent was added
    lastwell = fin[lastindex].well

    if(fin[lastindex].reag!=op.reag):
        cost=1
    else:
        # check on all 3 remaining reagent types
        cost=0
        for i in range(0,len(w[0])):
            if(w[lastwell][i]!=fin[lastindex].reag):
                if not (added[lastwell][i]==0 or (w[lastwell][i]==w[op.well][i] and added[op.well][i]==1)):
                    cost=1

        # take into account pipette capacity
        if(cap!=None):
            #only need to do that if cost is ostensibly 0 and the number of operations is less than the capacity
            if((cost==0) and (len(fin)>=cap)):
                # check if the next dose of vector doesn't fit into the pipette due to capacity limitations
                # to do that, see how many doses of current reagent have been delivered
                backforcap = 0
                while (backforcap<len(fin)):
                    backforcap += 1
                    if(fin[-backforcap].reag!=op.reag):
                        backforcap -= 1
                        break

                # if capacity of the current pipette tip with this reagent is exceeded, will have to change tip
                if(backforcap%cap==0):
                    cost=1

    return cost


# ---------------------------INPUT CONVERSION------------------------------------
"""
Input can be stored as:
1) a 2d-list w, 
2) list of all operations to do ops,
3) subsets of operations grouped by reagent (relevant for TSP method and reorderings)

Conversion between these data types is often necessary
"""

#read subsets from w
def w_to_subsets(w,subsets):
    for i in range(0, len(w)):
        for j in range(0, 4):
            match = False
            for k in range(0, len(subsets)):
                if (subsets[k].reag == w[i][j]):
                    match = True
                    break
            if (match):
                subsets[k].nuwell(i)
            else:
                subsets.append(Ss(w[i][j], i))

#read w from subsets
def subsets_to_w(subsets,w):
    #determine the number of wells
    maxwell=0
    for i in range(0, len(subsets)):
        locmax=max(subsets[i].wells)
        if (locmax>maxwell):
            maxwell=locmax
    maxwell+=1

    #initialise w
    emptyreags=['','','','']
    for i in range(0,maxwell):
        w.append(emptyreags.copy())

    #fill w
    for i in range(0, len(subsets)):
        for well in subsets[i].wells:
            w[well][ADDRESS[subsets[i].reag[0]]]=subsets[i].reag

#read a sequence of operations from w
def subsets_to_ops(subsets,ops):
    for i in range(0, len(subsets)):
        for j in range(0, len(subsets[i].wells)):
            ops.append(Oper(subsets[i].reag, subsets[i].wells[j]))

def ops_to_subsets(ops, subsets):
    for op in ops:
        ispresent = False
        for subset in subsets:
            if (op.reag == subset.reag):
                subset.nuwell(op.well)
                ispresent = True
                break
        if not ispresent:
            subsets.append(Ss(op.reag,op.well))


# -------------------------------JSON READER-------------------------------
# pass an empty w or subsets if you want them filled; pass None if not
def jsonreader(filename, w, subsets):
    # load file for reading
    jsonfile = open(filename, "r")
    input = json.load(jsonfile)

    ss = []  # initialise subsets
    dic = {'constructs': {}, 'reagents': {}}  # preset the dictionary
    reagclass = {'promoter': 'p', 'rbs': 'r', 'cds': 'c',
                 'terminator': 't'}  # preset indices of reagents to be recorded in the subsets list
    reagnum = {'p': 0, 'r': 0, 'c': 0, 't': 0}  # preset the number of reagents of a given class
    wellno = 0  # preset well counter

    # read the input
    for construct in input.items():
        # match construct to well number
        dic['constructs'][wellno] = construct[0]

        # get reagents
        parts = construct[1]['parts']
        for part in parts:
            if (part == 'backbone'):  # backbone does not count!
                continue
            # determine reagent name
            reagname = parts[part]['name']

            # check if such name is already in the dictionary
            match = False
            for prevname in dic['reagents'].values():
                if (reagname == prevname):
                    match = True
                    break

            if (match):  # if yes and we're dealing with subsets, just add the new operation to the subset it belongs to
                for ss_it in ss:
                    if (dic['reagents'][ss_it.reag] == reagname):
                        ss_it.nuwell(wellno)
            else:  # if no, update the dictionary
                nuentry = reagclass[part] + str(reagnum[reagclass[part]])  # determine which entry to put
                dic['reagents'][nuentry] = reagname  # put the entry into dictionary
                ss.append(Ss(nuentry, wellno))
                reagnum[reagclass[part]] += 1  # update number of reagents of this class

        wellno += 1

    if (subsets != None):  # record subsets
        subsets = ss
    if (w != None):  # record the well array
        subsets_to_w(ss, w)

    # return dictionary that allows to decode the input information from the outputs
    return dic


