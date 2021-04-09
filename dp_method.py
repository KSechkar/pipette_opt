# DP-BASED METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.1.0, 23.2.21

"""
Having received the array of wells and parts w,
for all additions to be performed records are created (different record for every position in the sequence of operations to do).
For every record, the best (yielding the lowest number of tip changes) previous operation is determined.
From these recordings, and optimised sequence of operations is reconstructed.

The resultant sequence of operations is stored in a class Oper list fin.
"""

import time

# import functions from own files
from input_generator import wgenerator
from auxil import *
from statespace_methods import getops


# ------------------------------CLASS DEFINITIONS----------------------------------
# recording fromat used in the DP algorithm; operation, best pre, best cost of operation sequence up to this addition
class DPrecord:
    # initialisation
    def __init__(self, op, w, pos):
        self.op=op # which operation is performed
        self.added=np.zeros((len(w), len(w[0])), dtype=bool) # state of wells at the moment (0 for unadded reagent, 1 for added); currently just initialised with zeros
        self.bestcost=-1 # cost of best operation sequence leading up to this operation; currently just initialised as -1
        self.previndex=-1 # index of the record for best prior operation; currenlty just initialised as -1
        self.changed=False # indicates if this operation needs a tip change
        self.pos=pos # position of this operation in the sequence

    # for printing the part type, destination well and which operation it is in the sequence
    def __str__(self):
        strRep = self.op.part + ' -> w' + str(self.op.well) + ' | operation no.' + str(self.pos)
        return strRep

# ----------------- SOLVER FUNCTION -------------------
def dp_method(w,fin,reord,caps):
    # PART 1: initial preparations

    # PART 1.1: get array of records for every operation on all possible positions
    dprecs = []
    getdprecs(w, dprecs, reord)

    # PART 2: get the sequence of operations

    # PART 2.1: deal with the records for the first operation in sequnce
    for rec in dprecs[0]:
        rec.bestcost = 1 # this is the first operation, so just 1 tip used
        rec.changed = True  # this is the first operation, so a new tip is needed
        rec.added[rec.op.well][rec.op.part[0]] = True # record that the operation in question has been made

    # PART 2.2: deal with all other records
    # consider the second operation, then the third, etc.
    for pos in range(1, len(dprecs)):
        #print(str(pos) + ' of ' + str(len(dprecs)) + ' operations')  # uncomment if need to track the progress

        # consider all records in this position
        for rec in dprecs[pos]:
            # find costs of making this operation the next after all possible prior operations
            dpcosts = []
            for maybeprev in dprecs[pos - 1]:
                dpcosts.append(getdpcost(rec, maybeprev, dprecs, w, caps))

            rec.bestcost = min(dpcosts) # find the one with the best cost
            rec.previndex = dpcosts.index(rec.bestcost) # record it as the previous operation
            prevrec = dprecs[pos - 1][rec.previndex] # get the record of the best-cost previous operation

            # copy the status of wells from the determined previous record and update it
            rec.added = prevrec.added.copy()
            rec.added[rec.op.well][rec.op.part[0]] = True

            # if needed, record that the tip must be changed here
            if (prevrec.bestcost < rec.bestcost):
                rec.changed = True


    # PART 3: get past

    # PART 3.1: find the record for the last operation with the lowest total cost
    frecs = dprecs[-1]
    finrec = min(frecs, key=lambda frecs: frecs.bestcost)

    # PART 3.2: write this operation in fin
    fin.append(finrec.op)
    fin[-1].changed = finrec.changed # indicate whether the tip has to be changed

    # PART 3.3: write all the prior operations into fin
    for pos in reversed(range(1, len(dprecs))):
        finrec = dprecs[finrec.pos - 1][finrec.previndex] # find previous operation
        fin.insert(0, finrec.op) # write it into fin
        fin[0].changed = finrec.changed # indicate whether the tip has to be changed



# ----------------- AUXILIARY FUNCTIONS -------------------
# get list of type DPOper operations from w
def getdprecs(w,dprecs,reord):
    # PART 1: get the list of all operations to be done
    ops=[]
    getops(w,ops,reord) # reord specifies if a reordering has to be applied (see auxil.py)


    # PART 2: create records for every operation on all possible positions in the sequence
    for i in range(0,len(ops)):
        dprecs.append([])
        for j in range(0,len(ops)):
            dprecs[i].append(DPrecord(ops[j], w, i))



# get cost of having the given operation after a potential previous one (1 if tip must be changed, 0 if not)
def getdpcost(rec,maybeprev,dprecs,w,caps):
    # if the operation has already been made, the resultant sequence is impossible
    if(maybeprev.added[rec.op.well][rec.op.part[0]]==True):
        return np.inf

    # if the part being added is different, need to change the tip
    if(rec.op.part!=maybeprev.op.part):
        extracost=1 # this is the cost of performing the given operation after the previous
    else:
        extracost=0

        # check if the last well had any parts the next well does not
        for i in range(0, len(w[0])):
            if (w[rec.op.well][i] != rec.op.part):
                if not (maybeprev.added[maybeprev.op.well][i] == 0 or (w[rec.op.well][i] == w[maybeprev.op.well][i] and maybeprev.added[rec.op.well][i] == 1)):
                    extracost=1

        # take into account pipette capacity, IF working with a capacitated problem
        if(caps!=None):
            # only need to do that if cost is ostensibly 0 and the number of operations is less than the capacity
            if(extracost==0 and caps[rec.op.part]<=rec.pos):
                extracost=1 # assume there is no extra space in the tip for another aliquot

                # check if it is actually so; if there is space, tip does not have to be changed
                checkrec=maybeprev
                for i in range(0,caps[rec.op.part]-1): # go back for as long as the capacity considerations are relevant
                    if(checkrec.changed):
                        extracost=0
                        break
                    else:
                        checkrec=dprecs[checkrec.pos-1][checkrec.previndex] # go to the record before

    # return the cost of the route up to the previous record plus the cost of performing the operation in question
    return extracost+maybeprev.bestcost


# ----------- MAIN FUNCTION (TESTING ONLY) ------------
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
    # get capacitites
    caps=capacities(reqvols,10,1.0)
    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # Call the solver. Specify the heuristic reordering used by changing reord
    # CHANGE THE 'REORD' ARGUMENT TO APPLY AN HEURSTIC REORDERING BEFORE SOLVING...
    dp_method(w, fin, reord=None, caps=caps)

    # display the solution
    dispoper(fin)

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')

    # calculate cost based on tip change indicators of the operations
    print('The total number of pipette tips used is (from resultant list) ' + str(route_cost(fin)))
    print('The total number of pipette tips used is (independent calculation) ' + str(independent_cost(fin, w, caps)[0]))

# main call
if __name__ == "__main__":
    main()