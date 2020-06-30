#TSP-BASED METHOD OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
#By Kirill Sechkar
#v0.0.4, 30.6.20

#The project makes use of the 'tspy' package

import numpy as np
import time
from tspy import TSP #TSP solver package
from tspy.solvers.utils import get_cost

#-------------------------------CLASS DEFINITIONS-------------------------------
#Each reagent matched with a subest of wells it is added to
class Ss:
    def __init__(self, reag, wellno): #initialisation
        self.reag=reag
        self.wells=[wellno]
        
    def nuwell(self, wellno): #record new well in the subset
        self.wells.append(wellno)
        
    def __str__(self): #for printing the subset's reagent type and wells out
        strRep=self.reag+'|'
        for i in range(0, len(self.wells)):
            strRep=strRep+' '+str(self.wells[i])
        return strRep
    
#Final output format is an array of Operations - will be the same for ALL methods
class Oper:
    def __init__(self, reag, well):
        self.reag=reag
        self.well=well
        
    def __str__(self): #for printing the subset's reagent type and wells out
        strRep=self.reag+' -> w'+str(self.well)
        return strRep
 
#-------------------------------INPUT-------------------------------
#Will be replaced by a test example generator or manual input reading function

#this is the example given to me in the main pipette_opt file
w=[['p1', 'r2', 'c4', 't2'],
   ['p2', 'r2', 'c1', 't2'],
   ['p1', 'r3', 'c2', 't1'],
   ['p2', 'r3', 'c1', 't1']]


#-------------------------------MAIN-------------------------------
def main():
    subsets=[] #array of all subsets (class Ss variables)
    fin=[] #final array where the operations are to be recorded
    tipchanges=0 #counts the total number of tip changes
    reagdic={} #dictionary that matches actual reagent names with p1, r2, c0, etc.
    
    #initialise the matrix of distances, i.e. our graph of wells
    D=np.zeros((len(w),len(w)))
    
    #Get the subsets:
    #option 1: read a .json file [comment to deselect]
    #jsonreader('level_zero_constructs.json',subsets,reagdic)
    #option 2: use a pre-set 2D list [comment to deselect]
    convert(w,subsets)
    
    tipchanges=len(subsets) #anyhow, we have to change the tip between the different reagents and we have a tip at first
    
    #print subsets and D (TEST ONLY)
    disp(subsets, D)
    
    #reorder the subsets. currently: in random order
    randreorder(subsets)
    
    #implement the algorithm
    for i in range(0,len(subsets)):
        tipchanges=singlesub(subsets[i],D,fin,tipchanges)
        
    dispoper(fin)  
    print('The total number of pipette tips used is '+str(tipchanges))
    
#-------------------------------FUNCTIONS-------------------------------  
def convert(w, subsets):  
    #create all subsets
    for i in range(0,len(w)):
        for j in range(0,4):
            match=False
            for k in range(0,len(subsets)):
                if(subsets[k].reag==w[i][j]):
                    match=True
                    break
            if(match):
                subsets[k].nuwell(i)
            else:
                subsets.append(Ss(w[i][j],i))
                
def disp(subsets, D):
    for i in range(0,len(subsets)):
        print(subsets[i])
    print(D)    

def dispoper(fin):
    for i in range(0,len(fin)):
        print(fin[i])

def jsonreader(filename,subsets,reagdic):
    jsonfile = open(filename,"r")
    if(jsonfile.mode=="r"):
        jsoncontent=jsonfile.readlines() #read the file into an array of text line strings
        
        #preset well indicies and indices of reagents to be recorded in the subsets list
        well=-1
        reagnum={'p': 0, 'r': 0, 'c': 0, 't': 0}
        
        jsonline=0
        while(jsonline<len(jsoncontent)):
            #if the new construct description starts, we assign it a new well number
            if(jsoncontent[jsonline][3:6]=='ss_'):
                well+=1
            
            print(jsonline)    
            #determine which reagent class we're to deal with    
            if(jsoncontent[jsonline][7:15]=='promoter'):
                reagclass='p'
            elif(jsoncontent[jsonline][7:10]=='rbs'):
                reagclass='r'
            elif(jsoncontent[jsonline][7:10]=='cds'):
                reagclass='c'
            elif(jsoncontent[jsonline][7:17]=='terminator'):
                reagclass='t'
            elif(jsoncontent[jsonline][7:15]=='backbone'): #if it's the backbone, skip
                jsonline+=3
            
            #determine reagent name
            if(jsoncontent[jsonline][9:13]=='name'):
                reagname=''
                for jsonletter in jsoncontent[jsonline][17:]:
                    if(jsonletter!='"'):
                        reagname+=jsonletter
                    else:
                        break
                print(reagname)
                #make a corresponding addition to subsets
                match=False
                for reagdic_it in reagdic:
                    if(reagname==reagdic[reagdic_it]):
                        match=True
                        break
                if(match):
                    for subsets_it in subsets:
                        if(subsets_it.reag==reagdic_it):
                            subsets_it.nuwell(well)
                else:
                    reagdic_newkey=reagclass+str(reagnum[reagclass])
                    subsets.append(Ss(reagdic_newkey,well))
                    reagdic[reagdic_newkey]=reagname
                    reagnum[reagclass]+=1
                
            jsonline+=1 #increase the counter
            
#-------------------------------REORDERINGS-------------------------------
#Will try various reorderings, not only just random

#!!Question: is it best to a)reshuffle the original array ('shuffle'); 
#                          b)copy subd]sets into a new reshuffle array ('permutation'); 
#                          c)create an array of reshuffled subsets addresses?

def randreorder(subsets):  #randomly reshuffle using default seed
    np.random.shuffle(subsets)
    
def randtimereorder(subsets):  #randomly reshuffle using time as seed
    np.random.RandomState(seed=round(time.time())).shuffle(subsets)
    
#-------------------------------SOLVE TSP FOR ONE SUBSET-------------------------------
def singlesub(subset,D,fin,tipchanges):
    #PART 1: initial preparations
    #get length to avoid calling len too often
    sublen=len(subset.wells)
    
    #initialise the subset's matrix subD
    subD=np.zeros((sublen+1,sublen+1)) #vertex 0, all edges to and from it being zero, allows to use cyclic TSP soluction for our PATH problem
    
    #PART 2: select submatrix and update D as if problem for the subset is already solved
    for i_well in range(0,sublen):
        current_well=0
        for j_D in range(0,len(D)):
            if(j_D==subset.wells[current_well]):
                subD[i_well+1][current_well+1]=D[subset.wells[i_well]][j_D] #select the edges within the subset into the submatrix
                if(current_well<sublen-1):
                    current_well+=1
            else:
                D[subset.wells[i_well]][j_D]=1 #make the edge going from the subset into the rest of D equal to one (updating D)
 
    #PART 3: solve TSP for the subset
    tsp=TSP()
    tsp.read_mat(subD)
    
    from tspy.solvers import TwoOpt_solver
    two_opt = TwoOpt_solver(initial_tour='NN', iter_num=100)
    tour = tsp.get_approx_solution(two_opt)
    print(tour)
    
    #PART 4: record the operations into the final output, 'unwrapping' the cycle arround the added zero node to create a path
    #find the position of the zero node in the tour
    i=0
    while(tour[i]!=0):
        i+=1
    #record the part after the zero node
    i+=1
    while(i<len(tour)-1):
        fin.append(Oper(subset.reag,subset.wells[tour[i]-1]))
        i+=1
    #record the part 'before'the zero node
    i=0
    while(tour[i]!=0):
        fin.append(Oper(subset.reag,subset.wells[tour[i]-1]))
        i+=1
        
    #PART 5: return the adjusted number of pipette tip changes
    return tipchanges+get_cost(tour,tsp) #include the tour cost in the number of tip changes

#-------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()
