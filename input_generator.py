#GENERATOR OF TEST INPUTS
#By Kirill Sechkar
#v0.0.1, 30.6.20

import numpy as np
import time

#simple test input generator. Out of given numbers of reagents of each class, create random input mixtures
def wgenerator(howmany_inputs, howmany_p, howmany_r, howmany_c, howmany_t):
    #initialise w and onewell
    w=[]
    
    #randomly fill w
    for i in range(0, howmany_inputs):
        w.append(['p'+str(np.random.randint(0,howmany_p)),
                  'r'+str(np.random.randint(0,howmany_r)),
                  'c'+str(np.random.randint(0,howmany_c)),
                  't'+str(np.random.randint(0,howmany_t))])
    
    return w
 
#main function  is only needed to test out the generation functions       
def main():
    print(wgenerator(2,4,2,6,2))

if __name__ == "__main__":
        main()