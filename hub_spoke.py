# HUB-AND-SPOKE METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.1, 27.7.20

import numpy as np
import time
from itertools import combinations, product

from input_generator import wgenerator
from auxil import dispoper

ADDRESS = {'p': 0, 'r': 1, 'c': 2, 't': 3}

# -------------------------------CLASS DEFINITIONS-------------------------------
class Oper:
    def __init__(self, reag, well):
        self.reag = reag
        self.well = well

    def __str__(self):  # for printing the subset's reagent type and wells out
        strRep = self.reag + ' -> w' + str(self.well)
        return strRep

class Hubwell:
    def __init__(self):
       self.wellno=[]
       self.content=[]
       self.origin=[]

    def __init__(self,content,wellno,origins):
       self.wellno=wellno
       self.content=content
       self.origins=origins

class Hub:
    def __init__(self,reags):
        self.reags = reags
        self.rank = len(reags)

        self.next=[0,1,2,3]
        for r in reags:
            self.next.pop(self.next.index(ADDRESS[r[0]]))

        self.wells=[]

    def getwell(self,well):
        self.wells.append(well)

    def __str__(self):
        strRep = 'Reagents: ' + str(self.reags)

        strRep += '\nWells: '
        for i in range(0,len(self.wells)):
            strRep += '\n' + str(self.wells[i])
        return strRep


# -------------------------------INPUT-------------------------------
# Will be replaced by a test example generator or manual input reading function
# this is the example given to me in the main pipette_opt file
w = [['p1', 'r2', 'c4', 't2'],
     ['p1', 'r2', 'c5', 't2'],
     ['p1', 'r2', 'c2', 't2'],
     ['p1', 'r2', 'c1', 't2']]


# -------------------------------MAIN-------------------------------
def main():
    fin = []  # an Oper list of operations in the order which they should be performed

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to dfine the number of wells
    change 4 last arguments to define the size of p, r, c and t reagent sets"""
    w = wgenerator(96, 6, 6, 3, 4)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    hubspoke(w,fin)
    # PERFORMANCE EVALUATION: print the working time
    dispoper(fin)
    print('The total number of pipette tip changes is '  + str(costfromhubs()))
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')


# -------------------------------SOLVER-------------------------------
def hubspoke(w,fin):
    # get list of reagents (sortent and unsorted by type)
    reaglist,listall=countreags(w)

    # get unfilled hubs
    global hubs
    hubs={0:{},1:{},2:{},3:{}}
    for n in reaglist.keys():
        hubs[1].update({(i): Hub(tuple([i])) for i in reaglist[n]})
    for m,n in combinations(reaglist.keys(),2):
        hubs[2].update({(i,j): Hub((i,j)) for i,j in product(reaglist[m],reaglist[n])})
    for m,n,o in combinations(reaglist.keys(),3):
        hubs[3].update({(i,j,k): Hub((i,j,k)) for i,j,k in product(reaglist[m],reaglist[n],reaglist[o])})

    # assign each well to some level-1 hub - randomly
    for i in range(0,len(w)):
        # randreag=np.random.randint(0,4)
        randreag=0
        hubs[1][w[i][randreag]].getwell(Hubwell(w[i],i,[0]))

    rankup(1)
    rankup(2)
    purgerank(3)


    # TEST ONLY
    ranked = 0
    for h in hubs[2].values():
        if(len(h.wells)!=0):
            ranked += len(h.wells)
    print(ranked)

    makefin(w,fin)


# -------------------------------FUNCTIONS USED BY SOLVER-------------------------------
def countreags(w):
    reaglist={'p':[],'r':[],'c':[],'t':[]}
    listall=[]
    for i in range(0, len(w)):
        for j in range(0, 4):
            match = False
            for r in reaglist[w[i][j][0]]:
                if (r == w[i][j]):
                    match = True
                    break
            if not (match):
                reaglist[w[i][j][0]].append(w[i][j])
                listall.append(w[i][j])
    return reaglist,listall


def rankup(rank):
    # transfer wells to level-2 hubs
    for hub in hubs[rank].values():
        delwell = []
        for i in range(0, len(hub.wells)):
            for j in range(0, len(hub.next)):
                nextaddress = addr(hub.reags + tuple([hub.wells[i].content[hub.next[j]]]))
                if ((len(hubs[rank + 1][nextaddress].wells) != 0) or (j == (len(hub.next) - 1))):
                    nuwell = hub.wells[i]
                    nuwell.origins.append(hub.reags)
                    hubs[rank + 1][nextaddress].getwell(nuwell)
                    delwell.append(i)
                    break
        if (len(delwell) != 0):
            dell = len(delwell)
            for d in range(1, dell + 1):
                hub.wells.pop(delwell[dell - d])


def purgerank(rank):
    for hub in hubs[rank].values():
        if (len(hub.wells) == 1):
            if(rank>2):
                origin = hub.wells[0].origins[-1]
            else:
                origin = hub.wells[0].origins[-1][0]

            nuwell=hub.wells[0]
            nuwell.origins.pop()

            hubs[rank-1][origin].getwell(nuwell)
            hub.wells.pop()
    if (rank>2):
        purgerank(rank-1)
    return 0


def addr(reags):
    ordered=[]
    for type in 'prct':
        for r in reags:
            if (r[0]==type):
                ordered.append(r)
                break
    if (len(ordered)==1):
        return tuple([ordered[0]])
    elif (len(ordered)==2):
        return (ordered[0],ordered[1])
    else:
        return (ordered[0],ordered[1],ordered[2])


def makefin(w,fin):
    for hub in hubs[1].values():
        if(len(hub.wells)!=0):
            for hubwell in hub.wells:
                fin.append(Oper(hub.reags[0],hubwell.wellno))
            for hubwell in hub.wells:
                for n in hub.next:
                    fin.append(Oper(w[hubwell.wellno][n],hubwell.wellno))

    multiwell = 100
    for rank in range(2,4):
        for hub in hubs[rank].values():
            if (len(hub.wells)!=0):
                multireag = ''
                for r in hub.reags:
                    fin.append(Oper(r,multiwell))
                    multireag += r
                multiwell+=1

                for hubwell in hub.wells:
                    fin.append(Oper(multireag,hubwell.wellno))

                for hubwell in hub.wells:
                    for n in hub.next:
                        fin.append(Oper(w[hubwell.wellno][n], hubwell.wellno))


def costfromhubs():
    cost = 0
    for rank in range(1, 4):
        coeff = 4 - rank
        for hub in hubs[rank].values():
            if(len(hub.wells)!=0):
                cost = cost + rank + len(hub.wells) * coeff
    return cost


# -----------------------ADVANCED SOLVER------------------------

# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
