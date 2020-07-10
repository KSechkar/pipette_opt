#GENERATOR OF TEST INPUTS
#By Kirill Sechkar
#v0.1.0, 10.7.20

import numpy as np
import csv
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

#create a .csv file listing randomly generated inputs
def inputlist(hm_inputs,hm_wells,hm_p,hm_r,hm_c,hm_t):
    #create file name
    filename=str(hm_inputs)+'i_'+str(hm_wells)+'w_'+str(hm_p)+'p_'+str(hm_r)+'r_'+str(hm_c)+'c_'+str(hm_t)+'t.csv'

    with open(filename, mode="w+",newline='') as infile:
        infile_write=csv.writer(infile,delimiter=',')
        for i in range(0,hm_inputs):
            w=wgenerator(hm_wells,hm_p,hm_r,hm_c,hm_t)
            for j in range(0,len(w)):
                infile_write.writerow([w[j][0],w[j][1],w[j][2],w[j][3]])
            infile_write.writerow(['new input:'])

#main function  is only needed to test out the generation functions       
def main():
    print(wgenerator(2,4,2,6,2))

if __name__ == "__main__":
        main()