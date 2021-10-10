# TESTING AND COMPARING VARIOUS METHODS
# By Kirill Sechkar
# v0.1.0, 1.6.21

import csv
import time
import statistics as stats
import argparse
import sys, os

from ppopt.statespace import nns, greedy_tree
from ppopt.lp import lp_method
from ppopt.dp import dp_method
from ppopt.auxil import *


#----------------------------PERFORM TESTING----------------------------
def test(inputname,which,read,mini,maxi):
    # make output directories if not present already
    if not (os.path.exists('results')):
        os.mkdir('results')
    if not (os.path.exists('times')):
        os.mkdir('times')
    if not (os.path.exists('progress')):
        os.mkdir('progress')

    # create a label to the filename describing the arguments
    labelfile = (inputname.split('/')[-1]).split('.')[0]

    # make column labels
    columns = ['Number of wells']
    for i in range(mini, maxi + 1):
        columns.append(str(i))

    # initialise results array with row labels
    means = [['Nearest Neighbour'], ['NNs depth2'], ['Greedy'],
             ['Nearest Neighbour+sametogether'], ['NNs depth2+sametogether'], ['Greedy+sametogether'],
             ['LP'], ['LP+random'], ['LP+sametogether'], ['LP+greedy'],
             ['DP'], ['DP+random'], ['DP+sametogether'], ['DP+leastout']]

    # determine which algorithms to run, change results array accordingly
    if (which == 'all'):
        torun = range(0, len(means))
    elif (which == 'statespace'):
        for k in range(0, 8):
            means.pop(-1)
        torun = range(0, 6)
    elif (which == 'LP'):
        for k in range(0, 6):
            means.pop(0)
        for k in range(0, 4):
            means.pop(-1)
        torun = range(0, 4)
    elif (which == 'LPnr'):
        for k in range(0, 6):
            means.pop(0)
        for k in range(0, 6):
            means.pop(-1)
        torun = range(0, 2)
    elif (which == 'DP'):
        for k in range(0, 10):
            means.pop(0)
        torun = range(0, 4)
    elif (which == 'DPnr'):
        for k in range(0, 10):
            means.pop(0)
        for k in range(0, 2):
            means.pop(-1)
        torun = range(0, 2)
    elif (which == 'DPsl'):
        for k in range(0, 12):
            means.pop(0)
        torun = range(0, 2)
    elif (which == 'DPLPn'):
        for k in range(0, 6):
            means.pop(0)
        for k in range(0, 3):
            means.pop(1)
        for k in range(0, 3):
            means.pop(-1)
        torun = range(0, 2)
    elif (which == 'ssno'):
        for k in range(0, 11):
            means.pop(-1)
        torun = range(0, 3)
    elif (which == 'sssametogether'):
        for k in range(0, 3):
            means.pop(0)
        for k in range(0, 8):
            means.pop(-1)
        torun = range(0, 3)
    elif (which == 'DPn'):
        for k in range(0, 10):
            means.pop(0)
        for k in range(0, 3):
            means.pop(-1)
        torun = [0]
    elif (which == 'DPr'):
        for k in range(0, 11):
            means.pop(0)
        for k in range(0, 2):
            means.pop(-1)
        torun = [0]
    elif (which == 'DPs'):
        for k in range(0, 12):
            means.pop(0)
        for k in range(0, 1):
            means.pop(-1)
        torun = [0]
    elif (which == 'DPl'):
        for k in range(0, 13):
            means.pop(0)
        torun = [0]
    elif (which == 'ssr'):
        means = [['Nearest Neighbour+random'], ['NNs depth2+random'], ['Greedy+random']]
        torun = range(0, 3)
    elif (which == 'ssl'):
        means = [['Nearest Neighbour+leastout'], ['NNs depth2+leastout'], ['Greedy+leastout']]
        torun = range(0, 3)
    else:
        print('Error! Unspecified selection of algorithms')
        exit(1)

    # initialise other results array with selected algorithm names
    medians = deepcopy(means)
    stdevs = deepcopy(means)
    timemeans = deepcopy(means)
    timedevs = deepcopy(means)
    # open file with inputs
    with open(inputname, mode="r") as infile:
        infile_read = csv.reader(infile)

        for i in range(mini,maxi+1):
            # skip the read rows until the inputs of size i start
            for row in infile_read:
                if (row[0] == '{} constructs:'.format(i)):
                    break

            # allocate space for storing outcomes and runtimes
            all_sols = []
            all_times = []
            for j in range(0, len(means)):
                all_sols.append([])
                all_times.append([])

            # read and process the specified number of inputs of size i
            for j in range(0, read):
                # get an input
                w = []
                onewell = [(), (), (), ()]
                for row in infile_read:
                    if (row[0] == 'end of input'):
                        break
                    else:
                        for entry in range(0, len(row)):
                            onewell[entry] = tuple(map(int, row[entry].split(',')))
                        w.append(onewell.copy())

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
                    if (means[itr][0][0:2] == 'LP'):
                        if (len(means[itr][0]) == 2):
                            timer = time.time()
                            lp_method(w, fin, reord=None, caps=caps, maxtime=1)
                            timer = time.time() - timer
                        else:
                            timer = time.time()
                            lp_method(w, fin, means[itr][0][3:], caps=caps, maxtime=1)
                            timer = time.time() - timer
                    elif (means[itr][0][0:2] == 'DP'):
                        if (len(means[itr][0]) == 2):
                            timer = time.time()
                            dp_method(w, fin, reord=None, caps=caps)
                            timer = time.time() - timer
                        else:
                            timer = time.time()
                            dp_method(w, fin, means[itr][0][3:], caps=caps)
                            timer = time.time() - timer
                    else:
                        # define reorderingy
                        if (means[itr][0][-12:] == 'sametogether'):
                            reord = 'sametogether'
                        elif (means[itr][0][-6:] == 'random'):
                            reord = 'random'
                        elif (means[itr][0][-8:] == 'leastout'):
                            reord = 'leastout'
                        else:
                            reord = None

                        # get solution
                        if (means[itr][0][:7] == 'Nearest'):
                            timer = time.time()
                            nns(w, fin, 1, reord, caps)
                            timer = time.time() - timer
                        elif (means[itr][0][:3] == 'NNs'):
                            timer = time.time()
                            nns(w, fin, 2, reord, caps)
                            timer = time.time() - timer
                        elif (means[itr][0][:6] == 'Greedy'):
                            timer = time.time()
                            greedy_tree(w, fin, 'optimistic+cap', reord, caps)
                            timer = time.time() - timer

                    # get route cost and record
                    rc = route_cost(fin)
                    all_sols[itr].append(rc)
                    all_times[itr].append(timer)

            # get means/medians/standard devioations and record
            for itr in range(0, len(means)):
                means[itr].append(str(stats.mean(all_sols[itr])))
                medians[itr].append(str(stats.median(all_sols[itr])))
                stdevs[itr].append(str(stats.stdev(all_sols[itr])))
                timemeans[itr].append(str(stats.mean(all_times[itr])))
                timedevs[itr].append(str(stats.stdev(all_times[itr])))

            if not (os.path.exists('progress')):
                os.mkdir('progress')
            with open('progress/' + labelfile + '_log.txt', mode="w+") as progress:
                progress.write('Case for ' + str(i) + ' wells processed - ' + str(maxi - i) + ' to go')

    # record results in output files
    # open files
    outmeans = open('results/' + labelfile + '_means.csv', mode="w+", newline='')
    outmeans_w = csv.writer(outmeans, delimiter=',')
    outmedians = open('results/' + labelfile + '_medians.csv', mode="w+", newline='')
    outmedians_w = csv.writer(outmedians, delimiter=',')
    outstdevs = open('results/' + labelfile + '_stdevs.csv', mode="w+", newline='')
    outstdevs_w = csv.writer(outstdevs, delimiter=',')
    outtimemeans = open('times/' + labelfile + '_timemeans.csv', mode="w+", newline='')
    outtimemeans_w = csv.writer(outtimemeans, delimiter=',')
    outtimedevs = open('times/' + labelfile + '_timedevs.csv', mode="w+", newline='')
    outtimedevs_w = csv.writer(outtimedevs, delimiter=',')

    # put column labels
    outmeans_w.writerow(columns)
    outmedians_w.writerow(columns)
    outstdevs_w.writerow(columns)
    for itr in range(0, len(means)):
        outmeans_w.writerow(means[itr])
        outmedians_w.writerow(medians[itr])
        outstdevs_w.writerow(stdevs[itr])
        outtimemeans_w.writerow(timemeans[itr])
        outtimedevs_w.writerow(timedevs[itr])
    return


#---------------------------ARGUMENT PARSING----------------------------
def getArgs(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parse arguments.")
    parser.add_argument("-w", "--which", help="Which alroithms to test.")
    parser.add_argument("-r", "--read", type=int, help="How many inputs to read.")
    parser.add_argument("-min", "--mini", type=int, help="From which number of wells to start.")
    parser.add_argument("-max", "--maxi", type=int, help="At which number of wells to finish.")
    parser.add_argument("-n", "--name", type=str, help="Name of file containing the inputs.")
    arguments = parser.parse_args(args)
    return arguments


# ---------------------------------MAIN---------------------------------
# which specifies which algorithms are tested
# read is how many inputs we read form each file
# mini is the maximum number of wells we test
# maxi is the maximum number of wells we test
def main():
    # get arguments from command line
    arguments=getArgs()
    which=getattr(arguments,'which')
    read = getattr(arguments, 'read')
    mini = getattr(arguments, 'mini')
    maxi=getattr(arguments,'maxi')
    inputname='inputs/'+getattr(arguments,'name')

    test(inputname,which,read,mini,maxi)
    return


# -------------------------------MAIN CALL------------------------------
if __name__ == "__main__":
    main()