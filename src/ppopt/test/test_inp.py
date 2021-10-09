# GENERATOR OF TEST INPUTS - simply to conveniently create inputs for 2 to 96 wells in one go
# By Kirill Sechkar
# v0.1.0, 1.6.21

from input_generator import inputlist
import os


# main function  is only needed to test out the generation functions
def main():
    if not (os.path.exists('inputs')):
        os.mkdir('inputs')
    for i in range(2, 97):
        inputlist('inputs/',100, i, 6, 6, 3, 4)


if __name__ == "__main__":
    main()
