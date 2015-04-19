#!/usr/bin/python
# This program creates a set of maternal plasma reads
# Presumptions: sqrt(ff) ~ norm(.123,.035)
# sqrt(ff)**2 * 100 = ff

from argparse import ArgumentParser
import cProfile
import copy
import random

from numpy.random import normal as normdist

from config import *


#Constants
READS = 1000000
LENCHR = 53000000
COVERAGE = random.randint(1, 7)

parser = ArgumentParser()
parser.add_argument("-t",
    metavar="Type",
    type=str,
    choices=["22q11del", "22q11dup"],#, "22q13del", "complete", "longd"],
    default="22q11dup",
    help="Type of aneuploidy 22")
parser.add_argument("-p",
    metavar="Parent",
    type=str,
    choices=["p", "m"],#, "22q13del", "complete", "longd"],
    default="p",
    help="Maternal or Paternal")

def loadGenomes():
    one = open("%sreads/4-17-21-3-38" % OUTPUTPATH, "r")
    two = open("%sreads/4-17-21-5-35" % OUTPUTPATH, "r")
    return one,two

def loadArray(genomeFile):
    pass

def del22q11(to_del, no_del):
    print "Simulating deleted 22q11..."
    #add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q11' not in z[2], to_del)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    
    #add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
        
    print "Done."
    return z

def dup22q11(to_dup, no_dup):
    print "Simulating duplicated 22q11..."
    #get only reads on 22q11
    x = filter(lambda x: 'q11' in x[2], to_dup)
    y = []
    #duplicate the reads, take care of coverage, from 1st parent
    for _ in range(COVERAGE * 2 - 1):#do not extend an extra time - want exactly n*COVERAGE*2 reads
        y.extend(x)
    x.extend(y)

    #add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q11' not in z[2], to_dup)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    z.extend(x)
    
    #add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_dup)
        
    print "Done."
    return z

def complete(g):
    pass

def del22q13(g):
    pass

def longd(g):
    pass

def noTrisomy(g):
    pass

def main(ff, type, parent):
    m,p = loadGenomes()
    maxp = int(sum(1 for line in p) - ff * READS)
    p.seek(random.randint(0,maxp))
    
    # Discard partial lines
    p.readline()
    
    # Parternal DNA
    paternal = [p.readline() for _ in range(int(ff * READS)/2)]
    paternal = [tuple(l[0:-2].split(",")) for l in paternal]

    # Maternal DNA in fetus
    maternal = [m.readline() for _ in range(int(ff * READS)/2)]
    maternal = [tuple(l[0:-2].split(",")) for l in maternal]


    if type == "22q11del" and parent == "p":
        fetal = del22q11(paternal, maternal)
    elif type == "22q11del" and parent == "m":
        fetal = del22q11(maternal, paternal)
    elif type == "22q11dup" and parent == "p":
        fetal = dup22q11(paternal, maternal)
    elif type == "22q11dup" and parent == "m":
        fetal = dup22q11(maternal, paternal)
   
    #fill the rest of the plasma reads
    g = copy.deepcopy(fetal)
    # Maternal DNA
    while(len(g) < READS):
        g.append(tuple(m.readline()[0:-2].split(",")))
    # Cleanup
    #g = [tuple(l[0:-2].split(",")) for l in g]


if __name__ == "__main__":
    print "Coverage:", COVERAGE
    args = parser.parse_args()
    main(normdist(FFMEAN,FFSTD)**2, args.t, args.p)
    #cProfile.run('main(normdist(FFMEAN,FFSTD)**2 * 100)')

