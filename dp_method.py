import time
from input_generator import wgenerator
from auxil import *
from copy import deepcopy

class DPrecord:
    def __init__(self, op, w):
        self.op=op
        self.opsfromnow = [op]  # is True if this operation involves a tip change
        self.added=np.ones((len(w), len(w[0])), dtype=bool)
        self.bestcost=-1 # CHANGED property must be different!

    def __str__(self):
        strRep = self.op.part + ' -> w' + str(self.op.well) + ' | ' + str(len(self.opsfromnow)) + ' op remaining'
        return strRep

# ----------------- SOLVER FUNCTION -------------------
def dp(w,fin,caps):
    dprecs=[]
    getdprecs(w,dprecs)
    address = addrfromw(w)

    # Loop!
    for rec in dprecs[len(dprecs)-1]:
        rec.bestcost=0
        rec.added[rec.op.well][address[rec.op.part[0]]]=False

    for pos in reversed(range(0,len(dprecs)-1)):
        for rec in dprecs[pos]:
            dpcosts=[]
            for maybenext in dprecs[pos+1]:
                dpcosts.append(getdpcost(rec,maybenext,w,address,caps))

            rec.bestcost = min(dpcosts)
            nextrec=dprecs[pos+1][dpcosts.index(rec.bestcost)]

            rec.opsfromnow =rec.opsfromnow + deepcopy(nextrec.opsfromnow)
            if (nextrec.bestcost < rec.bestcost):
                rec.opsfromnow[1].changed = True

            rec.added=deepcopy(nextrec.added)
            rec.added[rec.op.well][address[rec.op.part[0]]]=False

            # TEST ONLY
            trues=0
        print(pos)

    frecs=dprecs[0]
    bestrec=min(frecs,key=lambda frecs: frecs.bestcost)
    for op in bestrec.opsfromnow:
        fin.append(op)
    fin[0].changed=True


# ----------------- AUXILIARY FUNCTIONS -------------------
# get list of type DPOper operations from w
def getdprecs(w,dprecs):
    for i in range(0,len(w)*len(w[0])):
        dprecs.append([])
        for well in range(0, len(w)):
            for part in range(0, len(w[well])):
                dprecs[i].append(DPrecord(Oper(w[well][part], well), w))

# get cost of putting the current operation before the next suggested operation
def getdpcost(rec,maybenext,w,address,caps):
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

        if(extracost==0 and caps[rec.op.part]<len(maybenext.opsfromnow)):
            extracost=1
            for i in range(0,caps[rec.op.part]):
                if(maybenext.opsfromnow[i].changed):
                    extracost=0
                    break

    return extracost+maybenext.bestcost



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

    w = [['p1', 'r2', 'c4', 't1'],
         ['p1', 'r2', 'c1', 't1'],
         ['p1', 'r3', 'c2', 't1'],
         ['p1', 'r3', 'c1', 't1']]

    # w = wgenerator(96, 6, 6, 3, 4)

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

    # Call the solver. Specify the heuristic reordering used by changing reord
    dp(w, fin, caps)

    dispoper(fin)

    # PERFORMACE EVALUATION: print the working time
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')

    fincost = 0
    global changes
    changes=[]
    for op in fin:
        changes.append(int(op.changed))
        fincost += int(op.changed)

    print('The total number of pipette tips used is (from resultant list) ' + str(fincost))
    print('The total number of pipette tips used is (independent calculation) ' + str(route_cost_with_w(fin, w, caps)))

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
    return cost


# main call
if __name__ == "__main__":
    main()