'''
Created on May 1, 2015
'''

import argparse
from config import *
from os.path import isfile


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
    print "Extracting data..."
    if not isfile("%s%s" % (OUTPUTPATH, input_file)):
        print("Extraction target not found: %s" % input_file)
        return False
    
    with open("%s%s" % (OUTPUTPATH, input_file), "r") as input:
        with open("%s/reads/%s" % (OUTPUTPATH, output_file), "w") as output:
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
                output.write("%s,%s,%s\n" % (position,read,getLoc(position)))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.i, args.o)