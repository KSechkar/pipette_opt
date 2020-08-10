# HUB-AND-SPOKE METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.3, 7.8.20

import time
from itertools import combinations, product

from input_generator import wgenerator
from auxil import *

# -------------------------------CLASS DEFINITIONS-------------------------------
# Final output format is an array of Operations - will be the same for ALL methods
# wells of number <100 are a stand-in for 'hub' wells
class Oper:
    def __init__(self, part, well):
        self.part = part
        self.well = well

    def __str__(self):  # for printing the subset's part type and wells out
        strRep = self.part + ' -> w' + str(self.well)
        return strRep

# Hubwell is how a well that will be served by a hub is listed within a Hub variable
class Hubwell:
    def __init__(self):
       self.wellno=[] # number of the well in w
       self.content=[] # which parts the well contains
       self.origin=[] # the hubs it was previously served by

    def __init__(self,content,wellno,origins):
       self.wellno=wellno
       self.content=content
       self.origins=origins

# Hub is a well where n parts are mixed ( 1 <= n <= 3, 'hub rank' is the value of n)
# Hubs of rank 1 are the actual part source wells, so they're treated quite differently
class Hub:
    def __init__(self,parts,pipinfo):
        self.parts = parts # the part mixed in the hub
        self.rank = len(parts) # hub rank

        self.next=[] # which type parts we add added to a served well next (e.g. '0' for 'p' if hub has 'r,c,t')
        for i in range(0,len(globw[0])):
            self.next.append(i)
        typad = addrfromw(globw)  # get part type addresses from w
        for r in parts:
            self.next.pop(self.next.index(typad[r[0]]))

        self.wells=[] # list of wells
        self.sets={} # wells get promoted to a higher-ranked hub in sets if they can be transferred to the same hub
                     # here, dictionary key is the extra part in the new hub

        # CAPACITY DETERMINATION
        # find volume to be delivered to each spoke from this hub (assuming perfect mixing)
        delivol = 0
        for r in parts:
            delivol += pipinfo['reqvols'][r]

        self.cap = capac(pipinfo['pipcap'],delivol,pipinfo['airgap'])

    # add a new well served by the hub
    def getwell(self, well):
        self.wells.append(well) # add well

        # add well to all promotion sets it must belong to
        for n in self.next:
            part = well.content[n]
            nomatch = True
            for s in self.sets:
                if (s == part):
                    self.sets[s].append(len(self.wells) - 1)
                    nomatch = False
                    break
            if (nomatch):
                self.sets[part] = [(len(self.wells) - 1)]

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
        strRep = 'parts: ' + str(self.parts)

        strRep += '\nWells: '
        for i in range(0,len(self.wells)):
            strRep += '\n' + str(self.wells[i])
        return strRep


# -------------------------------INPUT-------------------------------
# Will be replaced by a test example generator or manual input reading function
w = [['p1', 'r2', 'c4','t1'],
     ['p2', 'r2', 'c1','t1'],
     ['p1', 'r3', 'c2','t2'],
     ['p2', 'r3', 'c1', 't1']]


# -------------------------------MAIN (TESTING ONLY)-------------------------------
def main():
    fin = []  # an Oper list of operations in the order which they should be performed

    """randomly generate w [comment to keep the hand-written example]
    change 1st argument to dfine the number of wells
    change 4 last arguments to define the size of p, r, c and t part sets"""
    w = wgenerator(96, 6, 6, 3, 4)

    # PERFORMACE EVALUATION: start the timer
    time1 = time.time()

    # define capacity-related information (not in testing, reqvols will be differently defined)
    pipinfo={'pipcap': 10, 'airgap': 1, 'reqvols':{}}

    # generate required volumes (for testing)
    ss = []
    w_to_subsets(w, ss)
    for s in ss:
        if (s.part[0] == 'p'):
            pipinfo['reqvols'][s.part] = 1.09
        elif (s.part[0] == 'r'):
            pipinfo['reqvols'][s.part] = 0.33
        elif (s.part[0] == 'c'):
            pipinfo['reqvols'][s.part] = 0.36
        else:
            pipinfo['reqvols'][s.part] = 0.75

    # call solver
    hubspoke(w,fin,pipinfo)

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

    # get list of parts (sortent and unsorted by type)
    partlist, listall = countparts()

    # determine the highest possible rank boundary (might change)
    hirank=len(partlist)+1

    if (hirank > 4): # only hubs of rank up to 3 are supported, so cut off any higher ranks
        hirank = 4

    # get unfilled hubs
    global hubs
    hubs={}
    for i in range(1,hirank): #initalise
        hubs[i]={}

    if(hirank>1): # create hubs of rank 1 if possible
        for n in partlist.keys():
            hubs[1].update({(i): Hub(tuple([i]),pipinfo) for i in partlist[n]})
    if(hirank>2): # create hubs of rank 2 if possible
        for m,n in combinations(partlist.keys(),2):
            hubs[2].update({(i,j): Hub((i,j),pipinfo) for i,j in product(partlist[m],partlist[n])})
    if(hirank>3): # create hubs of rank 3 if possible
        for m,n,o in combinations(partlist.keys(),3):
            hubs[3].update({(i,j,k): Hub((i,j,k),pipinfo) for i,j,k in product(partlist[m],partlist[n],partlist[o])})



    # assign each well to SOME level-1 hub
    for i in range(0,len(globw)):
        randpart=0
        hubs[1][globw[i][randpart]].getwell(Hubwell(globw[i],i,[0]))

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
# create a list of all parts used
def countparts():
    partlist = {} # list of all parts, grouped by type
    for t in globw[0]:
        partlist[t[0]]=[]

    listall=[]
    for i in range(0, len(globw)):
        for j in range(0, len(globw[0])):
            match = False
            for r in partlist[globw[i][j][0]]:
                if (r == globw[i][j]):
                    match = True
                    break
            if not (match):
                partlist[globw[i][j][0]].append(globw[i][j])
                listall.append(globw[i][j])
    return partlist,listall


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
                nexthublink = sorthublink(hub.parts + tuple([s]))
                nexthublen = len(hubs[rank + 1][nexthublink].wells)
                hublen = len(hub.wells)
                setlen = len(hub.sets[s])
                nexthubcap=hubs[rank + 1][nexthublink].cap
                tipsnow = hubtips(hublen,hub.cap,rank) + hubtips(nexthublen,nexthubcap,rank+1)
                tipsthen = hubtips(hublen - setlen,hub.cap,rank) + hubtips(nexthublen + setlen,nexthubcap,rank+1)
                thisimprovement = tipsnow - tipsthen

                # if this is the new most-improved player, record that
                if (thisimprovement > improvement):
                    improvement = thisimprovement
                    bestset = s
            # stop when no improvement can happen any more
            if (improvement == 0):
                break

            #promote best set
            nexthublink = sorthublink(hub.parts + tuple([bestset]))
            for iw in hub.sets[bestset]:
                nuwell = hub.wells[iw]
                nuwell.origins.append(hub.parts)
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


# link to a hub has a strict ordering of part types as outlined by addresses typad
# this function makes unordered liks ordered
def sorthublink(parts):
    ordered=[]
    typad=addrfromw(globw)
    for type in sorted(typad, key=typad.get):
        for r in parts:
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
                fin.append(Oper(hub.parts[0],hubwell.wellno))
            for hubwell in hub.wells:
                for n in hub.next:
                    fin.append(Oper(globw[hubwell.wellno][n],hubwell.wellno))

    multiwell = 100
    for rank in range(2,len(hubs)+1):
        for hub in hubs[rank].values():
            if (len(hub.wells)!=0):
                multipart = ''
                for r in hub.parts:
                    fin.append(Oper(r,multiwell))
                    multipart += r
                multiwell+=1

                for hubwell in hub.wells:
                    fin.append(Oper(multipart,hubwell.wellno))

                for hubwell in hub.wells:
                    for n in hub.next:
                        fin.append(Oper(globw[hubwell.wellno][n], hubwell.wellno))


# calculate costs based on hub output
def costfromhubs():
    cost = 0
    for rank in range(1, len(hubs)+1):
        for hub in hubs[rank].values():
            cost += hubtips(len(hub.wells), hub.cap, rank)
    return cost


def hubtips(hublen, hubcap, rank):
    if (hublen == 0): # no wells => no tips at all needed
        return 0
    elif(hubcap==0): # if there is no way to deliver even one dose in a pipette at all
        return np.inf
    else:
        cost = rank * np.ceil(hublen / hubcap)  # number of tips needed to make the hub mixture
        cost += np.ceil(hublen / hubcap) - 1  # no. of tips to deliver mixture. -1 as we reuse last tip from before
        coeff = 4 - rank
        cost += coeff * hublen  # number of tips to deliver remaining parts to spoke wells
        return cost


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
