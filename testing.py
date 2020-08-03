# TESTING AND COMPARING VARIOUS METHODS
# By Kirill Sechkar
# v0.0.3, 8.7.20

import csv
import time
import statistics as stats

from tsp_method import tsp_method
from statespace_methods import iddfs, greedy_tree
from input_generator import wgenerator, inputlist
from auxil import route_cost_with_w, capac, commoncapac


# --------------------------------MAIN---------------------------------
def main():
    MAXI = 97  # maximum number of wells we test+1
    MINI = 2  # maximum number of wells we test+1
    READ = 100  # how many inputs we read form each file

    # make column labels
    columns = ['Number of wells']
    for i in range(MINI, MAXI):
        columns.append(str(i))

    # initialise results arrays with row labels
    means = [['TSP-random'], ['TSP'], ['TSP-sametogether'],
               ['TSP-nearest neighbour'], ['TSP-iddfs depth 2'], ['TSP-leastout'],
               ['TSP-greedy'],
               ['Nearest Neighbour'], ['iddfs depth 2'], ['Greedy'],
               ['Nearest Neighbour+sametogether'], ['iddfs depth 2+sametogether'],
               ['Greedy+sametogether']]
    medians=[['TSP-random'], ['TSP'], ['TSP-sametogether'],
               ['TSP-nearest neighbour'], ['TSP-iddfs depth 2'], ['TSP-leastout'],
               ['TSP-greedy'],
               ['Nearest Neighbour'], ['iddfs depth 2'], ['Greedy'],
               ['Nearest Neighbour+sametogether'], ['iddfs depth 2+sametogether'],
               ['Greedy+sametogether']]
    stdevs=[['TSP-random'], ['TSP'], ['TSP-sametogether'],
               ['TSP-nearest neighbour'], ['TSP-iddfs depth 2'], ['TSP-leastout'],
               ['TSP-greedy'],
               ['Nearest Neighbour'], ['iddfs depth 2'], ['Greedy'],
               ['Nearest Neighbour+sametogether'], ['iddfs depth 2+sametogether'],
               ['Greedy+sametogether']]

    # pipette capacity
    # find least capacity of all to make it the common value
    # With 40fmol of each part and concentrations from 'Start-Stop Assembly Calculator' it'll be 5
    cap = commoncapac(pipcap=10,airgap=1,filename='input/doses.csv')

    #get results
    for i in range(MINI,MAXI):
        all_sols=[]
        for j in range(0,len(means)):
            all_sols.append([])

        # open file with inputs
        filename='input/100i_'+str(i)+'w_6p_6r_3c_4t.csv'
        with open(filename,mode="r") as infile:
            infile_read = csv.reader(infile)
            for j in range(0,READ):
                w = nextw(infile_read)
                for itr in range(0,len(means)):
                    fin = []
                    if(means[itr][0][0:3]=='TSP'):
                        if(len(means[itr][0])==3):
                            tsp_method(w,fin,reord=None,filename=None,cap=cap)
                        else:
                            tsp_method(w,fin,means[itr][0][4:],filename=None,cap=cap)
                    else:
                        #define reordering
                        if(means[itr][0][-12:]=='sametogether'):
                            reord='sametogether'
                        else:
                            reord=None

                        #get solution
                        if(means[itr][0][:7]=='Nearest'):
                            iddfs(w,fin,1,True,reord,cap)
                        elif(means[itr][0][:5]=='iddfs'):
                            iddfs(w,fin,2,True,reord,cap)
                        elif (means[itr][0][:6] == 'Greedy'):
                            greedy_tree(w, fin, 'optimistic+cap', reord,cap)

                    #get route cost and record
                    rc=route_cost_with_w(fin,w,cap)
                    all_sols[itr].append(rc)

        #get means/medians/standard devioations and record
        for itr in range(0,len(means)):
            means[itr].append(str(stats.mean(all_sols[itr])))
            medians[itr].append(str(stats.median(all_sols[itr])))
            stdevs[itr].append(str(stats.stdev(all_sols[itr])))

        with open('progress/log.txt',mode="w+") as progress:
            progress.write('Case for '+str(i)+' wells processed - '+str(96-i)+' to go')

    #record results in output files
    #open files
    outmeans = open('results/means.csv', mode="w+",newline='')
    outmeans_w = csv.writer(outmeans, delimiter=',')
    outmedians = open('results/medians.csv', mode="w+",newline='')
    outmedians_w = csv.writer(outmedians, delimiter=',')
    outstdevs = open('results/stdevs.csv', mode="w+",newline='')
    outstdevs_w = csv.writer(outstdevs, delimiter=',')

    #put column labels
    outmeans_w.writerow(columns)
    outmedians_w.writerow(columns)
    outstdevs_w.writerow(columns)
    for itr in range(0,len(means)):
        outmeans_w.writerow(means[itr])
        outmedians_w.writerow(medians[itr])
        outstdevs_w.writerow(stdevs[itr])


#------------------------------CREATING AND READING PREDEFINED INPUT LIST------------
#create a list
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

#read next input
def nextw(infile_read):
    w=[]
    for row in infile_read:
        if (row == ['end of input']):
            break
        else:
            w.append(row)

    return w


def runtest(filename,hm_inputs):
    #open file with inputs
    infile = open(filename, mode="r")
    infile_read = csv.reader(infile)
    #open file to write outputs in
    outfile= open('out_best_tsp.csv',mode="w+")
    outfile_write=csv.writer(outfile,delimiter=',')
    #initialise output list of working times
    allsolutions=[[]]
    reord='sametogether'
    # various_reorderings=['greedy']
    #read inputs one-by-one and get solutions
    for i in range(hm_inputs):
        time1=time.time()
        for j in range(0,len(allsolutions)):
            print(' '+str(j))
            fin=[]
            w = nextw(infile_read)

            """if(j==0):
                iddfs(w,fin,1,True,reord)
            elif(j==1):
                greedy_tree(w,fin,'optimistic',reord)
            else:"""
            cap = capac(10, 1.5, 1)
            tips=tsp_method(w,fin,reord,filename=None,cap=cap)

            #tsp_method(w,fin,various_reorderings[j],None)
            allsolutions[j].append(str(tips))
        time1=1000*(time.time()-time1)

    for s in allsolutions:
        outfile_write.writerow(s)


# -------------------------------MAIN CALL-------------------------------
if __name__ == "__main__":
    main()