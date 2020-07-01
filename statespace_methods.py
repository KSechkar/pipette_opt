#STATE SPACE-BASED METHODS (SHORTEST PATH/TREE SEARCH) OF SOLVING THE PIPETTE TIP CHANGES OPTIMISATION PROBLEM
#By Kirill Sechkar
#v0.0.1, 30.6.20

import numpy as np
import time

from input_generator import wgenerator

#-------------------------------CLASS DEFINITIONS-------------------------------
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
    ops=[] #an Oper list of operations to be performed
    fin=[] #an Oper list of operations in the order which they should be performed
    
    #w=wgenerator(96,6,6,3,4)
    
    time1=time.time()
    
    tips=iddfs(w,fin)
    dispoper(fin)
    print(tips)
    
    print(time.time()-time1)
    
#-------------------------------SOLVERS-------------------------------
def iddfs(w,fin):
    ops=[]
    getops(w,ops)
    np.random.shuffle(ops)
    tips=1
    all_operations=len(w)*len(w[0])
    
    fin.append(ops[0])
    ops.pop(0)
    
    while(len(fin)<all_operations):
        tips+=iddfs_oneiter(ops,fin)
        
    return tips
    
    
def iddfs_oneiter(ops,fin):
    lastindex=len(fin)-1 #number of performed operations before the last one
    
    #get information about the last well where a reagent was added
    lastwell=fin[lastindex].well
    lastwell_addedreags=[]
    for i in range(0,lastindex):
        if(fin[i].well==lastwell):
            lastwell_addedreags.append(fin[i])
    
    #determine the cost of each possible operation
    cost=np.ones(len(ops))
    for op in range(0,len(ops)):
        cost[op]=1
        if(ops[op].reag==fin[lastindex].reag):
            for i in range(0,lastindex):
                for j in range(0,len(lastwell_addedreags)):
                    if(fin[i].well==ops[op].well) and (fin[i].reag==lastwell_addedreags[j]):
                        lastwell_addedreags.pop(j)
                        break
                if(len(lastwell_addedreags)==0):
                    cost[op]=0
                    break
    
    #pick the least costly operation
    nextop=cost.tolist().index(min(cost))
    fin.append(ops[nextop])
    ops.pop(nextop)
    return min(cost)
                
                
            
    
#-------------------------------AUXILLIARY FUNCTIONS-------------------------------  
#get a list of all operations from w
def getops(w,ops):
    for well in range(0,len(w)):
        for reagent in range(0,len(w[well])):
            ops.append(Oper(w[well][reagent],well))

#display an Oper list
def dispoper(fin):
    for i in range(0,len(fin)):
        print(fin[i])    
        
#-------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()