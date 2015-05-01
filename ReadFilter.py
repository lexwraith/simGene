'''
Created on May 1, 2015
'''

import argparse
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

def main(input_file, output_file):
    with open(input_file, "r") as input:
        with open(output_file, "w") as output:
            for line in input:
                if "chr22" in line:
                    output.write(line)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.i, args.o)