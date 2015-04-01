#!/usr/bin/python
'''
Created on Mar 27, 2015
'''

from subprocess import call
from os.path import isfile
from multiprocessing import Process
import os
from datetime import datetime as dt
from time import sleep

#Constants
MAXNUMPROCS = 10
path = {}
path['REF'] = 'REF.fa'
path['ALT'] = 'ALT.fa'

def removeFiles(target):
    if isfile("data/%s.aln" % target):
        os.remove("data/%s.aln" % target)
    if isfile("data/%s.fq" % target):
        os.remove("data/%s.fq" % target)
    if isfile("data/%s.sam" % target):
        os.remove("data/%s.sam" % target)

def extractReads(target):
    print "Extracting data..."
    if not isfile("data/" + target + "_errFree.sam"):
        print("Extraction target not found: %s" % target)
        return False

    with open("data/" + target + "_errFree.sam") as input_file:
        with open("data/reads/" + target, "a") as output_file:
            for line in input_file:
                #skip header files
                if line[0] == "@":
                    continue
                line_components = line.split()
                if len(line_components) < 4:
                    continue
                position = line_components[3]
                read = line_components[9]
                output_file.write(position + "    " + read + "\n")
    os.remove("data/%s_errFree.sam" % target)
    return True

def main(seq, cov):
    # Sets filename to Month-Day-Hour-Minute-Seconds name
    filename = (str(x) for x in dt.now().timetuple()[1:6])
    filename = reduce(lambda x,y: x + "-" +  y, filename)

    # For some reason can't have multiple reads on same file
    indcopy = filename + ".cp"
    print("Copying %s to %s" % (path[seq], indcopy))
    call("cp data/%s data/%s" % (path[seq], indcopy), shell=True)

    print "Generating reads..."
    toCall = "art_illumina -ef -sam -i data/%s -l 25 -ss HS25 -f %s -o data/%s" % (indcopy, cov, filename)
    print(toCall)
    return_code = call(toCall, shell=True)
    os.remove("data/%s" % indcopy)

    if return_code != 0:
        print("Error occurred generating reads for %s" % filename)
        print(return_code)
        return

    print "Reads generated successfully"
    removeFiles(filename)
    result = extractReads(filename)

    if not result:
        print "Error occurred extracting the reads"
        return

    print "Reads have been extracted."

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("s", 
        metavar="Sequence", 
        type=str, 
        choices=["REF","ALT"], 
        help="REFerence or ALTernate sequence.")
parser.add_argument("c",
        metavar="Coverage",
        type=int,
        default=1,
        help="Read coverage over any position.")

parser.add_argument("-n", 
        metavar="Iterations", 
        type=int,
        default=1,
        help="Number of runs.")

if __name__ == '__main__':
    args = parser.parse_args()
    arg = "(%s, %s)" % (args.s, args.c)
    p = [Process(target=main, args=(args.s,args.c)) for _ in range(args.n)]
    for i in p:
        i.start()
        sleep(1)
