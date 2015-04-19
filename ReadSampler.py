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
    choices=["none", "22q11del", "22q11dup", "22q13del", "complete", "longd"],
    default="none",
    help="Type of aneuploidy 22")
parser.add_argument("-p",
    metavar="Parent",
    type=str,
    choices=["p", "m"],
    default="p",
    help="Maternal or Paternal")

def loadGenomes():
    one = open("%sreads/4-17-21-3-38" % OUTPUTPATH, "r")
    two = open("%sreads/4-17-21-5-35" % OUTPUTPATH, "r")
    return one,two

'''
Creates fetal data for 22q11 deletion.
Final fetal data will include all of no_del, and all of to_del besides anything in 22q11
'''
def del22q11(to_del, no_del):
    print "Simulating deleted 22q11..."
    # Add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q11' not in z[2], to_del)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
        
    print "Done."
    return z

'''
Creates fetal data for 22q11 duplication
Final fetal data includes all of no_dup, all of to_dup, with everything in to_dup that is 22q11 twice
'''
def dup22q11(to_dup, no_dup):
    print "Simulating duplicated 22q11..."
    # Get only reads on 22q11
    x = filter(lambda x: 'q11' in x[2], to_dup)
    y = []
    # Duplicate the reads, take care of coverage, from 1st parent
    for _ in range(COVERAGE * 2 - 1):# Do not extend an extra time - want exactly n*COVERAGE*2 reads
        y.extend(x)
    x.extend(y)

    # Add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q11' not in z[2], to_dup)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    z.extend(x)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_dup)
        
    print "Done."
    return z

def complete(to_dup, no_dup):
    print "Simulating complete duplicate"
    z = copy.deepcopy(to_dup)
    y = []
    # Duplicate the reads, take care of coverage, from 1st parent
    for _ in range(COVERAGE * 2 - 1):# Do not extend an extra time - want exactly n*COVERAGE*2 reads
        y.extend(to_dup)
    z.extend(y)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_dup)
        
    print "Done."
    return z

def del22q13(to_del, no_del):
    print "Simulating deleted 22q13..."
    # Add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q13' not in z[2], to_del)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
        
    print "Done."
    return z

def longd(to_del, no_del):
    print "Simulating complete deletion"
    z = []
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
    
    print "Done."
    return z

def noAneuploidy(maternal, paternal):
    print "Simulating no aneuploidy..."
    z = copy.deepcopy(maternal)
    z.extend(paternal)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    print "Done."
    return z

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


    # Get reads for the fetus alone
    if type == "22q11del" and parent == "p":
        fetal = del22q11(paternal, maternal)
    elif type == "22q11del" and parent == "m":
        fetal = del22q11(maternal, paternal)
    elif type == "22q11dup" and parent == "p":
        fetal = dup22q11(paternal, maternal)
    elif type == "22q11dup" and parent == "m":
        fetal = dup22q11(maternal, paternal)
    elif type == "22q13del" and parent == "p":
        fetal = del22q13(paternal, maternal)
    elif type == "22q13del" and parent == "m":
        fetal = del22q13(maternal, paternal)
    elif type == "complete" and parent == "p":
        fetal = complete(paternal, maternal)
    elif type == "complete" and parent == "m":
        fetal = complete(maternal, paternal)
    elif type == "longd" and parent == "p":
        fetal = longd(paternal, maternal)
    elif type == "longd" and parent == "m":
        fetal = longd(maternal, paternal)
    elif type == "none":
        fetal = noAneuploidy(maternal, paternal)
    # Fill the rest of the plasma reads
    g = copy.deepcopy(fetal)
    # Maternal DNA
    while(len(g) < READS):
        g.append(tuple(m.readline()[0:-2].split(",")))
    # Cleanup
    #g = [tuple(l[0:-2].split(",")) for l in g]
    print "Writing to output file " + type + "_" + parent
    with open(OUTPUTPATH + type + "_" + parent, "w") as f:
        for entry in g:
            f.write(entry + "\n")
    print "Done."


if __name__ == "__main__":
    print "Coverage:", COVERAGE
    args = parser.parse_args()
    main(normdist(FFMEAN,FFSTD)**2, args.t, args.p)
    #cProfile.run('main(normdist(FFMEAN,FFSTD)**2 * 100)')

