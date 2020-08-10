# AUXIL.PY
# Auxiliary functions (e.g. route cost calculator, operation printer, json file reader) used by both TSP- and State Space-based methods
# v0.1.2, 10.8.2020

import numpy as np
import json
import csv

# -----------------------class definitions-----------------------------
#operations
class Oper:
    def __init__(self, part, well):
        self.part = part
        self.well = well

    def __str__(self):  # for printing the subset's part type and wells out
        strRep = self.part + ' -> w' + str(self.well)
        return strRep

#operation subsets
class Ss:
    def __init__(self, part, wellno):  # initialisation
        self.part = part
        self.wells = [wellno]

    def nuwell(self, wellno):  # record new well in the subset
        self.wells.append(wellno)

    def __str__(self):  # for printing the subset's part type and wells out
        strRep = self.part + '|'
        for i in range(0, len(self.wells)):
            strRep = strRep + ' ' + str(self.wells[i])
        return strRep


# -----------------------RESULTS PRINTING-----------------------------
def dispoper(fin):
    for i in range(0, len(fin)):
        print(fin[i])


# -----------------------CAPACITY CALCULATOR-----------------------------
# need to clairfy if the last collection has an air gap after it
# currently assumed that NO
def capac(pipcap, dose, airgap):
    doseandgap = dose + airgap
    cap = int(pipcap / doseandgap)  # get how many collections followed by gap will fit
    if (dose <= (pipcap - cap * doseandgap)):  # see if a final collection without air gap would also fit
        cap += 1
    return cap

# determine how many doses of each part a pipette can hold (i.e. find capacities)
def capacities(reqvols, pipcap, airgap):
    cap={}
    for reqvol in reqvols.keys():
        cap[reqvol]=capac(pipcap,reqvols[reqvol],airgap)
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

    # find the number of the last well where a part was added
    lastwell = fin[lastindex].well
    # parts present in this well affect the cost and thus need to be determined
    lastwell_addedparts = []
    for i in range(0, lastindex):
        if (fin[i].well == lastwell):
            lastwell_addedparts.append(fin[i].part)

    # determine the cost of each possible operation
    cost = 1  # by default, it is 1 for any operation
    if (op.part == fin[lastindex].part):
        # this part of the program checks the conditions outlined in the 'Problem Representation' Document
        for i in range(0, lastindex):
            for j in range(0, len(lastwell_addedparts)):
                if (fin[i].well == op.well) and (fin[i].part == lastwell_addedparts[j]):
                    lastwell_addedparts.pop(j)
                    break
            if (len(lastwell_addedparts) == 0):
                break
        if (len(lastwell_addedparts) == 0):
            cost = 0
    return cost


# alterative way to determine operation cost by also knowing the the input well array
def route_cost_with_w(fin,w,caps):
    added = np.zeros((len(w), len(w[0]))) #tells which parts were added to which well

    # get addresses of part types in w
    address = addrfromw(w)

    #for the first operation in fin
    cost=1
    added[fin[0].well][address[fin[0].part[0]]]=1
    for i in range(1, len(fin)):
        cost += cost_func_with_w(fin[0:i], fin[i], w, added,caps)
        added[fin[i].well][address[fin[i].part[0]]] = 1
    return cost


def cost_func_with_w(fin,op,w,added,caps):
    # find the index of last for easier further referencing
    lastindex = len(fin) - 1
    if(lastindex<0): #if no previous operations have been performed, we obviously need to put on a tip
        return 1

    # find the num ber of the last well where a part was added
    lastwell = fin[lastindex].well

    if(fin[lastindex].part!=op.part):
        cost=1
    else:
        # check on all 3 remaining part types
        cost=0
        for i in range(0,len(w[0])):
            if(w[lastwell][i]!=fin[lastindex].part):
                if not (added[lastwell][i]==0 or (w[lastwell][i]==w[op.well][i] and added[op.well][i]==1)):
                    cost=1

        # take into account pipette capacity
        if(caps!=None):
            #only need to do that if cost is ostensibly 0 and the number of operations is less than the capacity
            if((cost==0) and (len(fin)>=caps[op.part])):
                # check if the next dose of vector doesn't fit into the pipette due to capacity limitations
                # to do that, see how many doses of current part have been delivered
                backforcap = 0
                while (backforcap<len(fin)):
                    backforcap += 1
                    if(fin[-backforcap].part!=op.part):
                        backforcap -= 1
                        break

                # if capacity of the current pipette tip with this part is exceeded, will have to change tip
                if(backforcap%caps[op.part]==0):
                    cost=1

    return cost


# ---------------------------INPUT CONVERSION------------------------------------
"""
Input can be stored as:
1) a 2d-list w, 
2) list of all operations to do ops,
3) subsets of operations grouped by part (relevant for TSP method and reorderings)

Conversion between these data types is often necessary
"""

#read subsets from w
def w_to_subsets(w,subsets):
    for i in range(0, len(w)):
        for j in range(0, len(w[0])):
            match = False
            for k in range(0, len(subsets)):
                if (subsets[k].part == w[i][j]):
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
    emptyparts=['','','','']
    for i in range(0,maxwell):
        w.append(emptyparts.copy())

    #fill w
    address=addrfromw(w)
    for i in range(0, len(subsets)):
        for well in subsets[i].wells:
            w[well][address[subsets[i].part[0]]]=subsets[i].part

#read a sequence of operations from w
def subsets_to_ops(subsets,ops):
    for i in range(0, len(subsets)):
        for j in range(0, len(subsets[i].wells)):
            ops.append(Oper(subsets[i].part, subsets[i].wells[j]))

def ops_to_subsets(ops, subsets):
    for op in ops:
        ispresent = False
        for subset in subsets:
            if (op.part == subset.part):
                subset.nuwell(op.well)
                ispresent = True
                break
        if not ispresent:
            subsets.append(Ss(op.part,op.well))


# -------------------------------JSON READER-------------------------------
# pass an empty w or subsets if you want them filled; pass None if not
# ignorelist contains names of part types to be ignored (like backbone)
def jsonreader(filename, w, subsets,ignorelist):
    # load file for reading
    jsonfile = open(filename, "r")
    input = json.load(jsonfile)

    ss = []  # initialise subsets
    dic = {'constructs': {}, 'parts': {}}  # preset the dictionary
    partclass = {}  # preset indices of parts to be recorded in the subsets list
    address = {} # addresses of part types in w
    partnum = {} # current number of each part type different species
    wellno = 0  # preset well counter
    sid='abcdefghijklmnopqrstuvwxyz' # one-letter ids standing in for part types
    i_sid = 0 # auxiliary
    i_addr=0

    # read the input
    for construct in input.items():
        # match construct to well number
        dic['constructs'][wellno] = construct[0]

        # get parts
        parts = construct[1]['parts']
        for part in parts.keys():
            # skip an ignored part type
            for ignore in ignorelist:
                if(part==ignore):
                    continue

            #  if this is the first entry, fill
            if(wellno==0):
                partclass[part]=sid[i_sid]
                partnum[sid[i_sid]]=0
                address[sid[i_sid]]=i_addr
                i_sid+=1
                i_addr+=1

            # determine partent name
            partname = parts[part]['name']

            # check if such name is already in the dictionary
            match = False
            for prevname in dic['parts'].values():
                if (partname == prevname):
                    match = True
                    break

            if (match):  # if yes and we're dealing with subsets, just add the new operation to the subset it belongs to
                for ss_it in ss:
                    if (dic['parts'][ss_it.part] == partname):
                        ss_it.nuwell(wellno)
            else:  # if no, update the dictionary
                nuentry = partclass[part] + str(partnum[partclass[part]])  # determine which entry to put
                dic['parts'][nuentry] = partname  # put the entry into dictionary
                ss.append(Ss(nuentry, wellno))
                partnum[partclass[part]] += 1  # update number of parts of this class

        wellno += 1 # proceeding to next well
    if (subsets != None):  # record subsets
        subsets = ss
    if (w != None):  # record the well array
        subsets_to_w(ss, w)

    # return dictionary that allows to decode the input information from the outputs
    return dic, address

def inalist(key,list):
    match = False
    for k in list:
        if(k==key):
            match = True
    return match

def addrfromw(w):
    address = {}
    if(len(w)!=0):
        if(len(w[0])!=0):
            for i in range(0,len(w[0])):
                address[w[0][i][0]] = i # assign the address to the first letter code
    return address

jsonreader('level_zero_constructs.json',w=None,subsets=[],ignorelist=['backbone'])
