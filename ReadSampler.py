#!/usr/bin/python
# This program creates a set of maternal plasma reads
# Presumptions: sqrt(ff) ~ norm(.35,.1)
# sqrt(ff)**2 * 100 = ff

from argparse import ArgumentParser
from numpy.random import normal as normdist
import random
import cProfile

#Constants
FFMEAN = .35
FFSTD = .1
READS = 1000000
LENCHR = 53000000

parser = ArgumentParser()

def loadGenomes():
    one = open("data/reads/4-1-14-32-42", "r")
    two = open("data/reads/4-1-14-34-48", "r")
    return one,two

def loadArray(genomeFile):
    pass

def main(ff):
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

if __name__ == "__main__":
    args = parser.parse_args()
    main(normdist(FFMEAN,FFSTD)**2)
    #cProfile.run('main(normdist(FFMEAN,FFSTD)**2 * 100)')

