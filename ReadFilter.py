'''
Created on May 1, 2015
'''

import argparse
from config import *

parser = argparse.ArgumentParser()
parser.add_argument("i", 
        metavar="Input", 
        type=str,  
        help="Input file to be filtered")
parser.add_argument("o",
        metavar="Output",
        type=str,
        default=1,
        help="Output file")

def getLoc(pos):
    pos = int(pos)
    for elem in LOC:
        if elem[0] < pos < elem[1]:
            return elem[2]

def main(input_file, output_file):
    with open(input_file, "r") as input:
        with open(output_file, "w") as output:
            for line in input:
                #skip header files
                if line[0] == "@":
                    continue
                line_components = line.split()
                if len(line_components) < 4:
                    continue
                chromosome = line_components[2]
                if chromosome != "chr22":
                    continue
                position = line_components[3]
                read = line_components[9]
                output_file.write("%s,%s,%s\n" % (position,read,getLoc(position)))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.i, args.o)