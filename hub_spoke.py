# HUB-AND-SPOKE METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
# By Kirill Sechkar
# v0.0.3, 7.8.20

import time
from itertools import combinations, product

from input_generator import wgenerator
from auxil import *


# ---------------------------------CLASS DEFINITIONS---------------------------------------
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
    def __init__(self,parts):
        self.parts = parts # the part mixed in the hub
        self.rank = len(parts) # hub rank
        self.wells = []  # list of wells assigned to the hub

        # Hub.next is addresses of type parts which are added to a spoke well after hub mxture
        self.next=[]
        for i in range(0,len(globw[0])):
            self.next.append(i)
        typad = addrfromw(globw)  # get part type addresses from w
        for r in parts:
            self.next.pop(self.next.index(typad[r[0]]))

        self.sets={} # wells get promoted to a higher-ranked hub in sets if they can be transferred to the same hub
                     # (here, a dictionary key is the extra part in the new hub compared to the current hub)

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

    # for printing the hubs out
    def __str__(self):
        strRep = 'parts: ' + str(self.parts)

        strRep += '\nWells: '
        for i in range(0,len(self.wells)):
            strRep += '\n' + str(self.wells[i])
        return strRep


# --------------------------------------SOLVER---------------------------------------------
# solver function, returns cost
def hubspoke(w, fin, pipinfo):
    # PART 1: initial preparations

    # PART 1.1: make global copy of w and caps (to avoid passing them many times between functions)
    global globw, globcaps
    globw = w
    globcaps = capacities(pipcap=pipinfo['pipcap'],airgap=pipinfo['airgap'],reqvols=pipinfo['reqvols'])

    # PART 1.2: get list of parts grouped by type
    partlist = countparts()

    # PART 1.3: determine the highest possible rank boundary (might change)
    hirank=len(partlist)+1
    if (hirank > 4): # only hubs of rank up to 3 are supported, so cut off any higher ranks
        hirank = 4

    # PART 1.4: get unfilled hubs
    global hubs
    hubs={}
    # initialise
    for i in range(1,hirank):
        hubs[i]={}

    if(hirank>1): # create hubs of rank 1 if possible
        for n in partlist.keys():
            hubs[1].update({(i): Hub(tuple([i])) for i in partlist[n]})
    if(hirank>2): # create hubs of rank 2 if possible
        for m,n in combinations(partlist.keys(),2):
            hubs[2].update({(i,j): Hub((i,j)) for i,j in product(partlist[m],partlist[n])})
    if(hirank>3): # create hubs of rank 3 if possible
        for m,n,o in combinations(partlist.keys(),3):
            hubs[3].update({(i,j,k): Hub((i,j,k)) for i,j,k in product(partlist[m],partlist[n],partlist[o])})


    # PART 2: assign wells to hubs

    # PART 2.1: assign each well to SOME level-1 hub
    for i in range(0,len(globw)):
        hubs[1][globw[i][0]].getwell(Hubwell(globw[i],i,[0]))

    # PART 2.2: promote wells to higher-rank hubs if this is beneficial
    for r in range(2,hirank):
        rankup(r-1)


    # PART 3: record operations
    makefin(fin)

    # return the cost
    return costfromhubs()


# -----------------------------------FUNCTIONS USED BY SOLVER------------------------------
# create a list of all parts, grouped by type
def countparts():
    # initialise
    partlist = {}
    for t in globw[0]:
        partlist[t[0]]=[]

    # fill
    for i in range(0, len(globw)):
        for j in range(0, len(globw[0])):
            # see if the current part is already recorded
            match = False
            for r in partlist[globw[i][j][0]]:
                if (r == globw[i][j]):
                    match = True
                    break

            # if not, record it
            if not (match):
                partlist[globw[i][j][0]].append(globw[i][j])

    return partlist


# promote wells to next-rank hubs
def rankup(rank):
    # do this for each hub of current rank
    for hub in hubs[rank].values():
        delwell = [] # the wells that'll be promoted will have to be deleted from current hub; create this deletion list

        # every time, promote the set that brings best results, until promotion does not improve the result anymore
        while (len(hub.sets) != 0):
            # find best set
            improvement = 0
            for s in hub.sets:
                # get by how much promoting the set would improve situation
                nexthublink = sorthublink(hub.parts + tuple([s]))
                nexthublen = len(hubs[rank + 1][nexthublink].wells)
                hublen = len(hub.wells)
                setlen = len(hub.sets[s])
                nexthubparts=hubs[rank + 1][nexthublink].parts
                tipsnow = hubtips(hublen,hub.parts,rank) + hubtips(nexthublen,nexthubparts,rank+1)
                tipsthen = hubtips(hublen - setlen,hub.parts,rank) + hubtips(nexthublen + setlen,nexthubparts,rank+1)
                thisimprovement = tipsnow - tipsthen

                # if this is the new most-improving set, record that
                if (thisimprovement > improvement):
                    improvement = thisimprovement
                    bestset = s

            # stop when no improvement can happen any more
            if (improvement == 0):
                break

            #promote the best-improving set
            nexthublink = sorthublink(hub.parts + tuple([bestset]))
            for iw in hub.sets[bestset]:
                nuwell = hub.wells[iw]
                nuwell.origins.append(hub.parts)
                hubs[rank + 1][nexthublink].getwell(nuwell)
                delwell.append(iw)

            # clear current set, remove promoted wells from remaining sets
            hub.setout(bestset)

        # clear hub of all wells that have been promoted
        delwell.sort(reverse=True)
        for d in delwell:
            hub.wells.pop(d)


# link to a hub has a strict ordering of part types as outlined by addresses typad
# this function makes unordered liks ordered
def sorthublink(parts):
    ordered = [] # intialise the ordered link
    typad = addrfromw(globw) # get addresses

    # make the link ordered
    for type in sorted(typad, key=typad.get):
        for r in parts:
            if (r[0]==type):
                ordered.append(r)
                break

    # return the ordered link
    if (len(ordered)==1):
        return tuple([ordered[0]])
    elif (len(ordered)==2):
        return (ordered[0],ordered[1])
    else:
        return (ordered[0],ordered[1],ordered[2])


# calculate costs based on hub output
def costfromhubs():
    cost = 0
    for rank in range(1, len(hubs)+1):
        for hub in hubs[rank].values():
            cost += hubtips(len(hub.wells), hub.parts, rank)
    return cost


def hubtips(hublen, hubparts, rank):
    if (hublen == 0): # no wells => no tips at all needed => simply return 0
        return 0
    else:
        cost = 0
        for i in range(0,rank):
            cost += np.ceil(hublen/globcaps[hubparts[i]]) # add costs of adding each part to hub mixture

        # mixture delivered by the last tip that was on => cost unchanged

        cost += (4 - rank) * hublen  # number of tips to deliver remaining parts to spoke wells
        return cost


# create final list of operations from the hubs
def makefin(fin):
    # as rank 1 hubs are just the original part wells, they are treated differently
    for hub in hubs[1].values():
        if(len(hub.wells)!=0):
            # transfer contents of rank 1 hub
            for hubwell in hub.wells:
                fin.append(Oper(hub.parts[0],hubwell.wellno))
            # add all remaining parts
            for hubwell in hub.wells:
                for n in hub.next:
                    fin.append(Oper(globw[hubwell.wellno][n],hubwell.wellno))

    # hubs of rank > 1
    multiwell = 100 # address of a hub well is denoted as wN, where N>100
    for rank in range(2,len(hubs)+1):
        for hub in hubs[rank].values():
            if (len(hub.wells)!=0):
                # preparing hub mixture
                multipart = '' # 'multipart' shows which parts are mixed in the hub
                for r in hub.parts:
                    fin.append(Oper(r,multiwell)) # transfer mixture components to the hub
                    multipart += r # add new component's name to multipart
                multiwell += 1 # next hub well has a number that's greater than last by 1

                # transfer contents of the hub
                for hubwell in hub.wells:
                    fin.append(Oper(multipart,hubwell.wellno))

                # add all remaining parts
                for hubwell in hub.wells:
                    for n in hub.next:
                        fin.append(Oper(globw[hubwell.wellno][n], hubwell.wellno))


# ---------------------------------MAIN (TESTING ONLY)------------------------------------
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
         ['p2', 'r2', 'c1', 't1'],
         ['p1', 'r3', 'c2', 't2'],
         ['p2', 'r3', 'c1', 't1']]

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


if __name__ == "__main__":
    main()
