#!/usr/bin/python
# This program creates a set of maternal plasma reads
# Presumptions: sqrt(ff) ~ norm(.123,.035)
# sqrt(ff)**2 * 100 = ff

from argparse import ArgumentParser
import cProfile
import copy
from operator import itemgetter
import random

import matplotlib.pyplot as plt

from numpy.random import normal as normdist
from scipy.stats import poisson

from config import *


#Constants
READS = 1000000
LENCHR = 53000000
BUCKET_SIZE = 1000
COVERAGE = random.randint(1, 7)

parser = ArgumentParser()
parser.add_argument("-t",
    metavar="Type",
    type=str,
    choices=["none", "22q11del", "22q11dup", "22q13del", "complete", "longd"],
    default="none",
    help="Type of aneuploidy 22")
parser.add_argument("-d", help="display read coverage in a graph", action='store_true')
parser.add_argument("-p",
    metavar="Parent",
    type=str,
    choices=["p", "m"],
    default="p",
    help="Maternal or Paternal")

def loadGenomes():
    one = open("%sreads/mother_filtered" % OUTPUTPATH, "r")
    two = open("%sreads/child_filtered" % OUTPUTPATH, "r")
    return one,two

'''
Creates fetal data for 22q11 deletion.
Final fetal data will include all of no_del, and all of to_del besides anything in 22q11
'''
def del22q11(to_del):
    print "Simulating deleted 22q11..."
    
    # Add in the reads not associated with 22q11, from 1st parent
    z_normal = filter(lambda z: 'q11' not in z[2], to_del)
    
    z = []
    z.extend(z_normal)
    
    
    '''
    # Add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q11' not in z[2], to_del)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
    '''
        
    print "Done."
    return z

'''
Creates fetal data for 22q11 duplication
Final fetal data includes all of no_dup, all of to_dup, with everything in to_dup that is 22q11 twice
'''
def dup22q11(to_dup):
    print "Simulating duplicated 22q11..."
    
    x = filter(lambda x: 'q11' in x[2], to_dup)
    y = filter(lambda y: 'q11' not in x[2], to_dup)
    
    z = []
    z.extend(x)
    z.extend(x)
    z.extend(y)
    
    '''
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
    '''
        
    print "Done."
    return z

def complete(to_dup):
    print "Simulating complete duplicate"
    y = copy.deepcopy(to_dup)
    z = []
    z.extend(y)
    z.extend(y)
    
    '''
    z = copy.deepcopy(to_dup)
    y = []
    # Duplicate the reads, take care of coverage, from 1st parent
    for _ in range(COVERAGE * 2 - 1):# Do not extend an extra time - want exactly n*COVERAGE*2 reads
        y.extend(to_dup)
    z.extend(y)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_dup)
        
    '''
    print "Done."
    return z

def del22q13(to_del):
    print "Simulating deleted 22q13..."
    
    z_normal = filter(lambda z: 'q11' not in z[2], to_del)
    
    z = []
    z.extend(z_normal)
    
    
    '''
    # Add in the reads not associated with 22q11, from 1st parent
    z = filter(lambda z: 'q13' not in z[2], to_del)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
        
    '''
    print "Done."
    return z

def longd(to_del):
    print "Simulating complete deletion"
    z = []
    '''
    # Add in all reads from the other parent
    for _ in range(COVERAGE):
        z.extend(no_del)
    '''
    print "Done."
    return z

def noAneuploidy(fetal):
    print "Simulating no aneuploidy..."
    print "Done."
    return fetal
    '''
    z = copy.deepcopy(maternal)
    z.extend(paternal)
    y = []
    for _ in range(COVERAGE - 1):
        y.extend(z)
    z.extend(y)
    
    print "Done."
    return z
    '''

'''
Given plasma data, creates an observed sequence. For each position, sees if the number of reads is greater or less than the expected amount
Put into buckets of length 10,000 each
'''
def getSequence(data):
    print len(data)
    print "Generating observed sequence..."
    #child = open("%sreads/child_filtered" % OUTPUTPATH, "r")
    num_reads = len(data)#len(list(child))
    expected_reads = num_reads*READ_LEN*BUCKET_SIZE/CHR_LEN
    
#    dist = poisson(expected_reads)
    low, high = poisson.interval(0.333, expected_reads)
    print low, high, expected_reads
    #pois = poisson(expected_reads)
    #low = pois*0.33
    #high = pois*0.667

#    low = expected_reads*0.8
#    high = expected_reads*1.2
    coverage = {}
    for read in data:
#        read_info = read.split(',')
        pos = int(read[0])
        read_len = len(read[1])
        for i in range(read_len):
            bucket = (pos + i)/BUCKET_SIZE
            coverage[bucket] = coverage.get(bucket, 0) + 1
    key_nums = {}
    for key in coverage:
        if coverage[key] < low:
            #print coverage[key]
            coverage[key] = "L"
            key_nums["L"] = key_nums.get("L", 0) + 1
        elif coverage[key] < high:
            coverage[key] = "N"
            key_nums["N"] = key_nums.get("N", 0) + 1
        else:
            #print coverage[key]
            key_nums["H"] = key_nums.get("H", 0) + 1
            coverage[key] = "H"
            
    print "to be low:", low
    print key_nums

    observed_seq = []
    #go through each dictionary entry in order
    #for i in range(max(coverage.keys())):
    #    observed_seq.append(coverage[i])
    for key in sorted(coverage):
        observed_seq.append(coverage[key])
    print "Done."
    return observed_seq

def displayCoverage(reads, type):
    print "Displaying coverage..."
    if type == "22q11del":
        desc = "Deletion on 22q11"
    elif type == "22q11dup":
        desc = "Trisomy on 22q11"
    elif type == "22q13del":
        desc = "Deletion on 22q13"
    elif type == "complete":
        desc = "Complete trisomy"
    elif type == "longd":
        desc = "Complete deletion"
    elif type == "none":
        desc = "No aneuploidy"
    
    values = [int(l[0]) for l in reads]
    
    plt.hist(values, bins=CHR_LEN/BUCKET_SIZE)
    plt.title(desc)
    plt.xlabel("Read position")
    plt.ylabel("Coverage")
    plt.show()

def main(ff, type, parent, display):
    m,f = loadGenomes()
    fetal = list(f)
    fetal = [tuple(l[0:-2].split(",")) for l in fetal]
    
    # Generate the aneuploidy for the entire fetus
    if type == "22q11del":
        fetal = del22q11(fetal)
    elif type == "22q11dup":
        fetal = dup22q11(fetal)
    elif type == "22q13del":
        fetal = del22q13(fetal)
    elif type == "complete":
        fetal = complete(fetal)
    elif type == "longd":
        fetal = longd(fetal)
    elif type == "none":
        fetal = noAneuploidy(fetal)
    
    # Get the fetus samples that will show up in the plasma
    fetal = random.sample(fetal, int(READS*ff))
    
    # Fill the rest of the plasma reads
    g = copy.deepcopy(fetal)
    # Maternal DNA
    while(len(g) < READS):
        g.append(tuple(m.readline()[0:-2].split(",")))
    
    #sort the data on read position
    sorted(g, key=itemgetter(0))
    
    # Get an observation sequence
    seq = getSequence(g)
    
    print "Writing to output file " + type + "_" + parent
    with open(OUTPUTPATH + type + "_" + parent, "w") as f:
        for entry in seq:
            f.write(entry + "\n")
    print "Done."
    
    if display:
        displayCoverage(g, type)
        
    '''
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
            f.write("%s,%s,%s\n" % (entry[0], entry[1], entry[2]))
    print "Done."
    
    '''


if __name__ == "__main__":
    #print "Coverage:", COVERAGE
    args = parser.parse_args()
    main(normdist(FFMEAN,FFSTD)**2, args.t, args.p, args.d)
    #cProfile.run('main(normdist(FFMEAN,FFSTD)**2 * 100)')

