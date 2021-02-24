# GENERATOR OF TEST INPUTS
# By Kirill Sechkar
# v0.1.1, 23.2.20

import numpy as np
import csv
import time
from input_generator import wgenerator, inputlist


# main function  is only needed to test out the generation functions
def main():
    for i in range(2, 97):
        inputlist('inputs/',100, i, 6, 6, 3, 4)


if __name__ == "__main__":
    main()
