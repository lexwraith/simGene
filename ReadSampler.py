#!/usr/bin/python
# This program creates a set of maternal plasma reads
# Presumptions: sqrt(ff) ~ norm(.123,.035)
# sqrt(ff)**2 * 100 = ff

from argparse import ArgumentParser
from numpy.random import normal as normdist
from config import *

import random
import cProfile

#Constants
READS = 1000000
LENCHR = 53000000

parser = ArgumentParser()
parser.add_argument("-t",
    metavar="Type",
    type=str,
    choices=["22q11del", "22q11dup", "22q13del", "complete", "longd"],
    default="complete",
    help="Type of aneuploidy 22")

def loadGenomes():
    one = open("%sreads/4-1-14-32-42" % OUTPUTPATH, "r")
    two = open("%sreads/4-1-14-34-48" % OUTPUTPATH, "r")
    return one,two

def loadArray(genomeFile):
    pass

def del22q11():
    pass

def dup22q11():
    pass

def complete():
    pass

def del22q13():
    pass

def longd():
    pass

def main(ff, type):
    m,p = loadGenomes()
    maxp = int(sum(1 for line in p) - ff * READS)
    p.seek(random.randint(0,maxp))
    
    # Discard partial lines
    p.readline()
    
    # Parternal DNA
    g = [p.readline() for _ in range(int(ff * READS)/2)]
    # Maternal DNA
    while(len(g) < READS):
        g.append(m.readline())
    # Cleanup
    g = [tuple(l[0:-2].split()) for l in g]
    
    if type == "22q11del":
        del22q11(g)
    elif type == "22q11dup":
        dup22q11(g)
    elif type == "22q13del":
        del22q13(g)
    elif type == "complete":
        complete(g)
    elif type == "longd":
        longd(g)

if __name__ == "__main__":
    args = parser.parse_args()
    main(normdist(FFMEAN,FFSTD)**2, args.t)
    #cProfile.run('main(normdist(FFMEAN,FFSTD)**2 * 100)')

