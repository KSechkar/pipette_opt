# HUB-AND-SPOKE METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.2, 28.7.20

import numpy as np
import time
from itertools import combinations, product

from input_generator import wgenerator
from auxil import dispoper, capac, addrfromw

# -------------------------------CLASS DEFINITIONS-------------------------------
# Final output format is an array of Operations - will be the same for ALL methods
# wells of number <100 are a stand-in for 'hub' wells
class Oper:
    def __init__(self, reag, well):
        self.reag = reag
        self.well = well

    def __str__(self):  # for printing the subset's reagent type and wells out
        strRep = self.reag + ' -> w' + str(self.well)
        return strRep

# Hubwell is how a well that will be served by a hub is listed within a Hub variable
class Hubwell:
    def __init__(self):
       self.wellno=[] # number of the well in w
       self.content=[] # which reagents the well contains
       self.origin=[] # the hubs it was previously served by

    def __init__(self,content,wellno,origins):
       self.wellno=wellno
       self.content=content
       self.origins=origins

# Hub is a well where n reagents are mixed ( 1 <= n <= 3, 'hub rank' is the value of n)
# Hubs of rank 1 are the actual reagent source wells, so they're treated quite differently
class Hub:
    def __init__(self,reags):
        self.reags = reags # the reagent mixed in the hub
        self.rank = len(reags) # hub rank


        self.next=[] # which type reagents we add added to a served well next (e.g. '0' for 'p' if hub has 'r,c,t')
        for i in range(0,len(globw[0])):
            self.next.append(i)
        typad = addrfromw(globw)  # get reagent type addresses from w
        for r in reags:
            self.next.pop(self.next.index(typad[r[0]]))

        self.wells=[] # list of wells
        self.sets={} # wells get promoted to a higher-ranked hub in sets if they can be transferred to the same hub
                     # here, dictionary key is the extra reagent in the new hub

    # add a new well served by the hub
    def getwell(self, well):
        self.wells.append(well) # add well

        # add well to all promotion sets it must belong to
        for n in self.next:
            reag = well.content[n]
            nomatch = True
            for s in self.sets:
                if (s == reag):
                    self.sets[s].append(len(self.wells) - 1)
                    nomatch = False
                    break
            if (nomatch):
                self.sets[reag] = [(len(self.wells) - 1)]

    # called when a given set has been promoted and now needs deletion
    # all the promoted wells are removed from remaining sets
    def setout(self,out):
        for iw in self.sets[out]:
            for n in self.next:
                for s in self.sets:
                    if ((s==self.wells[iw].content[n]) and (s!=out)):
                        set = self.sets[self.wells[iw].content[n]]
                        ind = set.index(iw)
                        set.pop(ind)
                        break
        self.sets.pop(out)


    def __str__(self):
        strRep = 'Reagents: ' + str(self.reags)

        strRep += '\nWells: '
        for i in range(0,len(self.wells)):
            strRep += '\n' + str(self.wells[i])
        return strRep


# -------------------------------INPUT-------------------------------
# Will be replaced by a test example generator or manual input reading function
w1 = [['p1', 'r2', 'c4','t1','f15'],
     ['p2', 'r2', 'c1','t1','f14'],
     ['p1', 'r3', 'c2','t2','f14'],
     ['p2', 'r3', 'c1', 't1', 'f14']]


# -------------------------------MAIN (TESTING ONLY)-------------------------------
def main():
    fin = []  # an Oper list of operations in the order which they should be performed

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to dfine the number of wells
    change 4 last arguments to define the size of p, r, c and t reagent sets"""
    w = wgenerator(96, 6, 6, 3, 4)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # define capacity-related information
    pipinfo={'pipcap': 10, 'onedose':1.1, 'airgap': 1}
    # call solver
    hubspoke(w1,fin,pipinfo)

    # PERFORMANCE EVALUATION: print the working time
    dispoper(fin)
    print('The total number of pipette tip changes is '  + str(costfromhubs()))
    print('The program took ' + str(1000 * (time.time() - time1)) + 'ms')


# -------------------------------SOLVER-------------------------------
# solver function, returns cost
def hubspoke(w, fin, pipinfo):
    # make global copy of w (to avoid passing it many times between functions)
    global globw
    globw = w

    # get list of reagents (sortent and unsorted by type)
    reaglist, listall = countreags()

    # determine the highest possible rank boundary (might change)
    hirank=len(reaglist)+1
    if (hirank > 4): # only hubs of rank up to 3 are supported, so cut off any higher ranks
        hirank = 4

    # get capaicty constraints for each rank,
    global caps
    caps = {}
    if (pipinfo != None):
        for i in range(1, hirank):
            caps[i] = capac(pipinfo['pipcap'], i * pipinfo['onedose'], pipinfo['airgap'])  # get capacity

            if (caps[i] == 0):  # zero capacity => is the ACTUAL boundary (higher-rank hubs can't deliver to spokes)
                hirank = i
                break
    else:
        for i in range(1, hirank):
            caps[i] = 96

    # get unfilled hubs
    global hubs
    hubs={}
    for i in range(1,hirank): #initalise
        hubs[i]={}

    if(hirank>1): # create hubs of rank 1 if possible
        for n in reaglist.keys():
            hubs[1].update({(i): Hub(tuple([i])) for i in reaglist[n]})
    if(hirank>2): # create hubs of rank 2 if possible
        for m,n in combinations(reaglist.keys(),2):
            hubs[2].update({(i,j): Hub((i,j)) for i,j in product(reaglist[m],reaglist[n])})
    if(hirank>3): # create hubs of rank 3 if possible
        for m,n,o in combinations(reaglist.keys(),3):
            hubs[3].update({(i,j,k): Hub((i,j,k)) for i,j,k in product(reaglist[m],reaglist[n],reaglist[o])})



    # assign each well to SOME level-1 hub
    for i in range(0,len(globw)):
        randreag=0
        hubs[1][globw[i][randreag]].getwell(Hubwell(globw[i],i,[0]))

    for r in range(2,hirank):
        rankup(r-1)

    """
    TEST ONLY
    ranked = 0
    for rank in range(1,4):
        for h in hubs[rank].values():
            if(len(h.wells)!=0):
                ranked += len(h.wells)
    print(ranked)
    """

    # record operations
    makefin(fin)

    return costfromhubs()


# -------------------------------FUNCTIONS USED BY SOLVER-------------------------------
# create a list of all reagents used
def countreags():
    reaglist = {} # list of all reagents, grouped by type
    for t in globw[0]:
        reaglist[t[0]]=[]

    listall=[]
    for i in range(0, len(globw)):
        for j in range(0, len(globw[0])):
            match = False
            for r in reaglist[globw[i][j][0]]:
                if (r == globw[i][j]):
                    match = True
                    break
            if not (match):
                reaglist[globw[i][j][0]].append(globw[i][j])
                listall.append(globw[i][j])
    return reaglist,listall


# promote wells to next-level hubs
def rankup(rank):
    # do this for each hub of current rank
    for hub in hubs[rank].values():
        delwell = [] # the well that'll be promoted will have to be deleted from current hub
        # every time, promote the set that brings best results, until no sets are left
        while (len(hub.sets) != 0):
            # find best set
            improvement = 0
            for s in hub.sets:
                # get by how much promoting the set would improve situation
                nexthublink = sorthublink(hub.reags + tuple([s]))
                nexthublen = len(hubs[rank + 1][nexthublink].wells)
                hublen = len(hub.wells)
                setlen = len(hub.sets[s])
                tipsnow = hubtips(hublen,rank) + hubtips(nexthublen,rank+1)
                tipsthen = hubtips(hublen - setlen,rank) + hubtips(nexthublen + setlen,rank+1)
                thisimprovement = tipsnow - tipsthen

                # if this is the new most-improved player, record that
                if (thisimprovement > improvement):
                    improvement = thisimprovement
                    bestset = s
            # stop when no improvement can happen any more
            if (improvement == 0):
                break

            #promote best set
            nexthublink = sorthublink(hub.reags + tuple([bestset]))
            for iw in hub.sets[bestset]:
                nuwell = hub.wells[iw]
                nuwell.origins.append(hub.reags)
                hubs[rank + 1][nexthublink].getwell(nuwell)
                delwell.append(iw)
            # clear current set, remove promoted wells from remaining sets
            hub.setout(bestset)

        # clear hub of all wells promoted
        delwell.sort(reverse=True)
        for d in delwell:
            hub.wells.pop(d)


# purge the highest rank of one-well hubs
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


# link to a hub has a strict ordering of reagent types as outlined by addresses typad
# this function makes unordered liks ordered
def sorthublink(reags):
    ordered=[]
    typad=addrfromw(globw)
    for type in sorted(typad, key=typad.get):
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


# create final list of operations using hub output
def makefin(fin):
    for hub in hubs[1].values():
        if(len(hub.wells)!=0):
            for hubwell in hub.wells:
                fin.append(Oper(hub.reags[0],hubwell.wellno))
            for hubwell in hub.wells:
                for n in hub.next:
                    fin.append(Oper(globw[hubwell.wellno][n],hubwell.wellno))

    multiwell = 100
    for rank in range(2,len(hubs)+1):
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
                        fin.append(Oper(globw[hubwell.wellno][n], hubwell.wellno))


# calculate costs based on hub output
def costfromhubs():
    cost = 0
    for rank in range(1, len(hubs)+1):
        for hub in hubs[rank].values():
            cost += hubtips(len(hub.wells), rank)
    return cost


def hubtips(hublen, rank):
    if (hublen == 0):
        return 0
    else:
        cost = rank * np.ceil(hublen / caps[1])  # number of tips needed to make the hub mixture
        cost += np.ceil(hublen / caps[rank]) - 1  # no. of tips to deliver mixture. -1 as we reuse last tip from before
        coeff = 4 - rank
        cost += coeff * hublen  # number of tips to deliver remaining reagents to spoke wells
        return cost


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
