# GENERATOR OF TEST INPUTS
# By Kirill Sechkar
# v1.0.0, 22.2.20

# Generate random Start-Stop Assembly inputs

import numpy as np
import csv
import time


# simple test input generator. Out of given numbers of parts of each class, create random input mixtures
def wgenerator(howmany_inputs, howmany_p, howmany_r, howmany_c, howmany_t):
    # initialise w and onewell
    w = []

    # randomly fill w
    for i in range(0, howmany_inputs):
        w.append([(0, np.random.randint(0, howmany_p)),
                  (1, np.random.randint(0, howmany_r)),
                  (2, np.random.randint(0, howmany_c)),
                  (3, np.random.randint(0, howmany_t))])

    return w


# create a .csv file listing randomly generated inputs
def inputlist(dir,hm_inputs, hm_wells, hm_p, hm_r, hm_c, hm_t):
    # create file name
    filename = str(hm_inputs) + 'i_' + str(hm_wells) + 'w_' + str(hm_p) + 'p_' + str(hm_r) + 'r_' + str(
        hm_c) + 'c_' + str(hm_t) + 't.csv'

    with open(dir+filename, mode="w+", newline='') as infile:
        infile_write = csv.writer(infile, delimiter=',')
        for i in range(0, hm_inputs):
            w = wgenerator(hm_wells, hm_p, hm_r, hm_c, hm_t)
            for j in range(0, len(w)):
                infile_write.writerow([str(w[j][0][0])+','+str(w[j][0][1]),
                                       str(w[j][1][0])+','+str(w[j][1][1]),
                                       str(w[j][2][0])+','+str(w[j][2][1]),
                                       str(w[j][3][0])+','+str(w[j][3][1])])
            infile_write.writerow(['end of input'])

def alldiff(dir,hm_wells,hm_types):
    filename = 'All_different_' + str(hm_wells) + '_wells.csv'
    with open(dir+filename, mode="w+", newline='') as infile:
        infile_write = csv.writer(infile, delimiter=',')
        for i in range(0, hm_wells):
            well_entry = []
            for j in range(0,hm_types):
                well_entry.append((j,i))
            infile_write.writerow(well_entry)
        infile_write.writerow(['end of input'])

def allsame(dir,hm_wells,hm_types):
    filename = 'All_same_' + str(hm_wells) + '_wells.csv'
    with open(dir+filename, mode="w+", newline='') as infile:
        infile_write = csv.writer(infile, delimiter=',')
        for i in range(0, hm_wells):
            well_entry = []
            for j in range(0,hm_types):
                well_entry.append((j,0))
            infile_write.writerow(well_entry)
        infile_write.writerow(['end of input'])

# main function  is only needed to test out the generation functions
def main():
    allsame('inputs/',96,4)


if __name__ == "__main__":
    main()
