# TESTING AND COMPARING VARIOUS METHODS
# By Kirill Sechkar
# v0.0.3, 8.7.20

import csv
import time
from copy import deepcopy
import statistics as stats

from tsp_method import tsp_method
from statespace_methods import nns, greedy_tree
from input_generator import wgenerator
from auxil import *


# --------------------------------MAIN---------------------------------
# which specifies which algorithms are tested
# read is how many inputs we read form each file
# mini is the maximum number of wells we test
# maxi is the maximum number of wells we test
def main(which,read,mini,maxi):
    # create a label to the filename describing the arguments
    labelfile = which + '_'+str(read)+'_' + str(mini) + '-' + str(maxi) + '_'

    # make column labels
    columns = ['Number of wells']
    for i in range(mini, maxi):
        columns.append(str(i))

    # initialise results array with row labels
    means = [['Nearest Neighbour'], ['NNs depth2'], ['Greedy'],
             ['Nearest Neighbour+sametogether'], ['NNs depth2+sametogether'], ['Greedy+sametogether'],
             ['LP'], ['LP+random'], ['LP+sametogether'], ['LP+greedy'],
             ['DP'], ['DP+random'], ['DP+sametogether']]

    # determine which algorithms to run, change results array accordingly
    if(which=='all'):
        torun=range(0,len(means))
    elif(which=='statespace'):
        for k in range(0,7):
            means.pop(-1)
        torun=range(0,6)
    elif(which=='statespace'):
        for k in range(0,6):
            means.pop(0)
        for k in range(0,3):
            means.pop(-1)
        torun=range(0,4)
    elif(which==''):
        for k in range(0,10):
            means.pop(0)
        torun=range(0,3)
    else:
        print('Error! Unspecified selection of algorithms')
        exit(1)

    # initialise other results array with selected algorithm names
    medians = deepcopy(means)
    stdevs = deepcopy(means)
    timemeans=deepcopy(means)
    timedevs=deepcopy(means)

    #get results
    for i in range(mini,maxi+1):
        all_sols=[]
        all_times=[]
        for j in range(0,len(means)):
            all_sols.append([])
            all_times.append([])

        # open file with inputs
        filename='inputs/100i_'+str(i)+'w_6p_6r_3c_4t.csv'
        with open(filename,mode="r") as infile:
            infile_read = csv.reader(infile)
            for j in range(0,read):
                w = nextw(infile_read)

                # generate required volumes (for testing). Values taken from a real instance of Start-Stop assembly
                ss = []
                w_to_subsets(w, ss)
                reqvols = {}
                for s in ss:
                    if (s.part[0] == 0):
                        reqvols[s.part] = 1.09
                    elif (s.part[0] == 1):
                        reqvols[s.part] = 0.33
                    elif (s.part[0] == 2):
                        reqvols[s.part] = 0.36
                    else:
                        reqvols[s.part] = 0.75

                # get capacitites
                caps = capacities(reqvols, 10, 1.0)

                for itr in torun:
                    fin = []
                    if(means[itr][0][0:2]=='LP'):
                        if(len(means[itr][0])==2):
                            timer=time.time()
                            tsp_method(w,fin,reord=None, caps=caps)
                            timer=time.time()-timer
                        else:
                            timer = time.time()
                            tsp_method(w,fin,means[itr][0][3:],caps=caps)
                            timer = time.time() - timer
                    elif(means[itr][0][0:2]=='DP'):
                        if (len(means[itr][0]) == 2):
                            timer = time.time()
                            tsp_method(w, fin, reord=None, caps=caps)
                            timer = time.time() - timer
                        else:
                            timer = time.time()
                            tsp_method(w, fin, means[itr][0][3:], caps=caps)
                            timer = time.time() - timer
                    else:
                        #define reordering
                        if(means[itr][0][-12:]=='sametogether'):
                            reord='sametogether'
                        else:
                            reord=None

                        # get solution
                        if(means[itr][0][:7]=='Nearest'):
                            timer = time.time()
                            nns(w,fin,1,reord,caps)
                            timer = time.time() - timer
                        elif(means[itr][0][:3]=='NNs'):
                            timer = time.time()
                            nns(w,fin,2,reord,caps)
                            timer = time.time() - timer
                        elif (means[itr][0][:6] == 'Greedy'):
                            timer = time.time()
                            greedy_tree(w, fin, 'optimistic+cap', reord,caps)
                            timer = time.time() - timer

                    # get route cost and record
                    rc=route_cost(fin,)
                    all_sols[itr].append(rc)
                    all_times[itr].append(timer)

        # get means/medians/standard devioations and record
        for itr in range(0,len(means)):
            means[itr].append(str(stats.mean(all_sols[itr])))
            medians[itr].append(str(stats.median(all_sols[itr])))
            stdevs[itr].append(str(stats.stdev(all_sols[itr])))
            timemeans[itr].append(str(stats.mean(all_times[itr])))
            timedevs[itr].append(str(stats.mean(all_times[itr])))

        with open('progress/'+labelfile+'log.txt',mode="w+") as progress:
            progress.write('Case for '+str(i)+' wells processed - '+str(maxi-i)+' to go')

    # record results in output files
    # open files
    outmeans = open('results/'+labelfile+'means.csv', mode="w+",newline='')
    outmeans_w = csv.writer(outmeans, delimiter=',')
    outmedians = open('results/'+labelfile+'medians.csv', mode="w+",newline='')
    outmedians_w = csv.writer(outmedians, delimiter=',')
    outstdevs = open('results/'+labelfile+'stdevs.csv', mode="w+",newline='')
    outstdevs_w = csv.writer(outstdevs, delimiter=',')
    outtimemeans= open('times/'+labelfile+'timemeans.csv', mode="w+",newline='')
    outtimemeans_w = csv.writer(outtimemeans, delimiter=',')
    outtimedevs = open('times/' + labelfile + 'timedevs.csv', mode="w+", newline='')
    outtimedevs_w = csv.writer(outtimedevs, delimiter=',')

    #put column labels
    outmeans_w.writerow(columns)
    outmedians_w.writerow(columns)
    outstdevs_w.writerow(columns)
    for itr in range(0,len(means)):
        outmeans_w.writerow(means[itr])
        outmedians_w.writerow(medians[itr])
        outstdevs_w.writerow(stdevs[itr])
        outtimemeans_w.writerow(timemeans[itr])
        outtimedevs_w.writerow(timedevs[itr])


#------------------------------CREATING AND READING PREDEFINED INPUT LIST------------
#read next input
def nextw(infile_read):
    w=[]
    onewell=[(),(),(),()]
    for row in infile_read:
        if (row == ['end of input']):
            break
        else:
            for entry in range(0,len(row)):
                onewell[entry] = tuple(map(int,row[entry].split(',')))
            w.append(onewell)

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
                nns(w,fin,1,True,reord)
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
    main('all',2,2,3)