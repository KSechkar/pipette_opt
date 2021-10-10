# GENERATOR OF TEST INPUTS - conveniently creates inputs for 2 to 96 wells in one go, similarly to what has been tested
# By Kirill Sechkar
# v0.1.0, 1.6.21

from input_generator import wgenerator
import os,csv


# main function - change the lines before calling test_inp to define which test inputs you require
def main():
    # directory where the inputs file will be saved
    dir='inputs/'
    # number of repeats for every input size
    howmany_inputs=100
    # numbers of possible promoters, RBSs, CDSs and terminators respectively
    howmany_p = 6
    howmany_r = 6
    howmany_c = 3
    howmany_t = 4

    if not (os.path.exists(dir[:-1])):
        os.mkdir(dir[:-1])
    test_inp(dir,howmany_inputs,howmany_p, howmany_r, howmany_c, howmany_t)
    return


def test_inp(dir, howmany_inputs,howmany_p, howmany_r, howmany_c, howmany_t):
    # name under which to save the test inputs
    filename = '{}Test_Inputs-{}i_{}p_{}r_{}c_{}t.csv'.format(dir,howmany_inputs,howmany_p, howmany_r,howmany_c, howmany_t)

    # generate and record a set of inputs for every input size
    with open(filename, 'w+', newline='') as infile:
        infile_write = csv.writer(infile, delimiter=',')
        for howmany_wells in range(2, 97):
            infile_write.writerow(['{} constructs:'.format(howmany_wells)])
            for i in range(0, howmany_inputs):
                w = wgenerator(howmany_wells, howmany_p, howmany_r, howmany_c, howmany_t)
                for j in range(0, len(w)):
                    infile_write.writerow([str(w[j][0][0]) + ',' + str(w[j][0][1]),
                                       str(w[j][1][0]) + ',' + str(w[j][1][1]),
                                       str(w[j][2][0]) + ',' + str(w[j][2][1]),
                                       str(w[j][3][0]) + ',' + str(w[j][3][1])])
                infile_write.writerow(['end of input'])
            infile_write.writerow([''])
    return
            

if __name__ == "__main__":
    main()
