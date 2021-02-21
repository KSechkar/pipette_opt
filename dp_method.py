import time
from copy import deepcopy

from input_generator import wgenerator
from auxil import *
from statespace_methods import getops

class DPrecord:
    def __init__(self, op, w, pos):
        self.op=op
        self.added=np.ones((len(w), len(w[0])), dtype=bool)
        self.bestcost=-1 # CHANGED property must be different!
        self.nextindex=-1
        self.nextchanged=False
        self.pos=pos

    def __str__(self):
        strRep = self.op.part + ' -> w' + str(self.op.well) + ' | operation no.' + str(self.pos)
        return strRep

# ----------------- SOLVER FUNCTION -------------------
def dp(w,fin,caps,reord):
    dprecs=[]
    getdprecs(w,dprecs,reord)
    address = addrfromw(w)

    # Loop!
    for rec in dprecs[len(dprecs)-1]:
        rec.bestcost=0
        rec.added[rec.op.well][address[rec.op.part[0]]]=False

    for pos in reversed(range(0,len(dprecs)-1)):
        for rec in dprecs[pos]:
            dpcosts=[]
            for maybenext in dprecs[pos+1]:
                dpcosts.append(getdpcost(rec,maybenext,dprecs,w,address,caps))

            rec.bestcost = min(dpcosts)
            rec.nextindex=dpcosts.index(rec.bestcost)

            nextrec = dprecs[pos + 1][rec.nextindex]
            if (nextrec.bestcost < rec.bestcost):
                rec.nextchanged = True

            rec.added=nextrec.added.copy()
            rec.added[rec.op.well][address[rec.op.part[0]]]=False

        print(pos)

    frecs=dprecs[0]
    finrec=min(frecs,key=lambda frecs: frecs.bestcost)
    fin.append(finrec.op)
    fin[0].changed=True
    for pos in range(1,len(dprecs)):
        nextchanged=finrec.nextchanged
        finrec=dprecs[finrec.pos+1][finrec.nextindex]
        fin.append(finrec.op)
        fin[-1].changed=nextchanged


# ----------------- AUXILIARY FUNCTIONS -------------------
# get list of type DPOper operations from w
def getdprecs(w,dprecs,reord):
    ops=[]
    getops(w,ops,reord)
    if(reord==None):
        for i in range(0,len(ops)):
            dprecs.append([])
            for j in range(0,len(ops)):
                dprecs[i].append(DPrecord(ops[j], w, i))
    else:
        for i in range(0,len(ops)):
            dprecs.append([])
            for j in reversed(range(0,len(ops))):
                dprecs[i].append(DPrecord(deepcopy(ops[j]), w, i))

# get cost of putting the current operation before the next suggested operation
def getdpcost(rec,maybenext,dprecs,w,address,caps):
    numops=len(w)*len(w[0])

    if(maybenext.added[rec.op.well][address[rec.op.part[0]]]==False):
        return np.inf

    if(rec.op.part!=maybenext.op.part):
        extracost=1
    else:
        extracost=0
        for i in range(0, len(w[0])):
            if (w[rec.op.well][i] != rec.op.part):
                if not (maybenext.added[rec.op.well][i] == 0 or (w[rec.op.well][i] == w[maybenext.op.well][i] and maybenext.added[maybenext.op.well][i] == 1)):
                    extracost=1

        if(extracost==0 and caps[rec.op.part]<numops-rec.pos):
            extracost=1
            checkrec=maybenext
            for i in range(1,caps[rec.op.part]):
                if(checkrec.nextchanged):
                    extracost=0
                    break
                else:
                    checkrec=dprecs[checkrec.pos+1][checkrec.nextindex]

    return extracost+maybenext.bestcost


def dpforward(w,fin,caps,reord):
    dprecs=[]
    getdprecs(w,dprecs,reord)
    address = addrfromw(w)

    for recs in dprecs:
        for rec in recs:
            rec.added=np.zeros((len(w),len(w[0])),dtype=bool)

    # Loop!
    for rec in dprecs[0]:
        rec.bestcost=0
        rec.added[rec.op.well][address[rec.op.part[0]]]=True
        rec.nextchanged=True

    for pos in range(1,len(dprecs)):
        for rec in dprecs[pos]:
            dpcosts=[]
            for maybeprev in dprecs[pos-1]:
                dpcosts.append(getdpcostforward(rec,maybeprev,dprecs,w,address,caps))

            rec.bestcost = min(dpcosts)
            rec.nextindex=dpcosts.index(rec.bestcost)

            prevrec = dprecs[pos - 1][rec.nextindex]
            if (prevrec.bestcost < rec.bestcost):
                rec.nextchanged=True

            rec.added=prevrec.added.copy()
            rec.added[rec.op.well][address[rec.op.part[0]]]=True
        print(pos)

    frecs=dprecs[-1]
    finrec=min(frecs,key=lambda frecs: frecs.bestcost)
    fin.append(finrec.op)
    fin[-1].changed=finrec.nextchanged
    for pos in reversed(range(1,len(dprecs))):
        finrec=dprecs[finrec.pos-1][finrec.nextindex]
        fin.insert(0,finrec.op)
        fin[0].changed=finrec.nextchanged

def getdpcostforward(rec,maybeprev,dprecs,w,address,caps):
    numops=len(w)*len(w[0])

    if(maybeprev.added[rec.op.well][address[rec.op.part[0]]]==True):
        return np.inf

    if(rec.op.part!=maybeprev.op.part):
        extracost=1
    else:
        extracost=0
        for i in range(0, len(w[0])):
            if (w[rec.op.well][i] != rec.op.part):
                if not (maybeprev.added[maybeprev.op.well][i] == 0 or (w[rec.op.well][i] == w[maybeprev.op.well][i] and maybeprev.added[rec.op.well][i] == 1)):
                    extracost=1

        if(extracost==0 and caps[rec.op.part]<=rec.pos):
            extracost=1
            checkrec=maybeprev
            for i in range(0,caps[rec.op.part]-1):
                if(checkrec.nextchanged):
                    extracost=0
                    break
                else:
                    checkrec=dprecs[checkrec.pos-1][checkrec.nextindex]

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

    w = [['p1', 'r2', 'c1', 't1'],
         ['p1', 'r2', 'c1', 't2'],
         ['p1', 'r2', 'c1', 't1'],
         ['p1', 'r2', 'c1', 't2']]

    w = wgenerator(96, 6, 6, 3, 4)

    # generate required volumes (for testing). Values taken from a real instance of Start-Stop assembly
    ss=[]
    w_to_subsets(w,ss)
    reqvols = {}
    for s in ss:
        if(s.part[0]=='p'):
            reqvols[s.part]=1.09
        elif(s.part[0]=='r'):
            reqvols[s.part]=0.33
        elif (s.part[0] == 'c'):
            reqvols[s.part] = 0.36
        else:
            reqvols[s.part] = 0.75
    # get capacitites
    caps=capacities(reqvols,10,1.0)
    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()
    #"""
    # Call the solver. Specify the heuristic reordering used by changing reord
    dp(w, fin, caps,reord=None)

    dispoper(fin)

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')
    
    fincost = 0
    for op in fin:
        fincost += int(op.changed)
    print('The total number of pipette tips used is (from resultant list) ' + str(fincost))
    route_cost_with_w(fin, w, caps)
    #"""
    fin2=[]
    dpforward(w, fin2, caps, reord=None)
    dispoper(fin2)
    ffincost = 0
    for op in fin2:
        ffincost += int(op.changed)
    print(ffincost)
    route_cost_with_w(fin2, w, caps)
    """
    vanillacost = 0
    for op in vanillafin:
        vanillacost += int(op.changed)
    print(vanillacost)
    """
    print(1)

def route_cost_with_w(fin,w,caps):
    # PART 1: initial preparations
    # get addresses of part types in w
    address = addrfromw(w)

    # make a copy, reset its tip change indicators
    cfin=deepcopy(fin)
    for i in range(0, len(cfin)):
        cfin[i].changed = False

    # create the array added (tells which parts were added to which well)
    added = np.zeros((len(w), len(w[0])),dtype=bool) # added[i][j]==True if part at address j in well i has been added



    # PART 2: get the cost

    # PART 2.1: the first operation in fin
    cost=1 # beginning the distribuiton => new tip taken
    cfin[0].changed = True # indicate the tip's been changed
    added[cfin[0].well][address[cfin[0].part[0]]] = 1 # indicate the part's been added

    # PART 2.2: all other operations
    for i in range(1, len(cfin)):
        one_cost = cost_func_with_w(cfin[0:i], cfin[i], w, added, caps) # get operation cost
        cost += one_cost # add operation cost

        added[cfin[i].well][address[cfin[i].part[0]]] = 1 # indicate the part's been added
        # if the tip's been changed (operation cost 1), indicate that
        if(one_cost==1):
            cfin[i].changed = True
        if (fin[i].changed != cfin[i].changed):
            print(str(cfin[i])+' || '+str(i))
    # PART 3: return the cost
    print('The total number of pipette tips used is (independent calculation) ' + str(cost))
    #dispoper(cfin)

# main call
if __name__ == "__main__":
    main()