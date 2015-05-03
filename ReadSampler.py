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
from _collections import defaultdict


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
    three = open("%sreads/father_filtered" % OUTPUTPATH, "r")
    return one,two, three

'''
Creates fetal data for 22q11 deletion.
Final fetal data will include all of the data called to the one parent, and only non-22q11 from the other
'''
def del22q11(to_del, parent):
    print "Simulating deleted 22q11..."
    
    # Add in the reads not associated with 22q11, from 1st parent
    z_normal = filter(lambda z: 'q11' not in z[2], to_del)
    # Add in the reads on 22q11, from the other parent
    z_q11 = filter(lambda z: 'q11' in z[2] and z[3] != parent and z[3] != "-", to_del)
        
    z = []
    z.extend(z_normal)
    z.extend(z_q11)
        
    print "Done."
    return z

'''
Creates fetal data for 22q11 duplication
Final data includes all non-22q11, 22q11 from one parent, and 2x 22q11 from the other
'''
def dup22q11(to_dup, parent):
    print "Simulating duplicated 22q11..."
    
    x = filter(lambda x: 'q11' in x[2] and x[3] != parent and x[3] != "-", to_dup)
    x_dup = filter(lambda x: 'q11' in x[2] and x[3] == parent or x[3] == "-", to_dup)
    y = filter(lambda y: 'q11' not in x[2], to_dup)
    
    z = []
    z.extend(x)
    z.extend(x_dup)
    z.extend(x_dup)
    z.extend(y)
        
    print "Done."
    return z

'''
Creates fetal data for a complete trisomy
Final data is all reads for one parent, and 2x reads for the other
'''
def complete(to_dup, parent):
    print "Simulating complete duplicate"
    x = filter(lambda x: x[3] != parent and x[3] != "-", to_dup)
    x_dup = filter(lambda x: x[3] == parent or x[3] == "-", to_dup)
    
    z = []
    z.extend(x)
    z.extend(x_dup)
    z.extend(x_dup)
    
    print "Done."
    return z

'''
Creates fetal data for a deleted 22q13
Final data is all reads from one parent, and all non-22q13 reads from the other
'''
def del22q13(to_del, parent):
    print "Simulating deleted 22q13..."
    
    # Add in the reads not associated with 22q11, from 1st parent
    z_normal = filter(lambda z: 'q13' not in z[2], to_del)
    # Add in the reads on 22q11, from the other parent
    z_q11 = filter(lambda z: 'q13' in z[2] and z[3] != parent and z[3] != "-", to_del)
        
    z = []
    z.extend(z_normal)
    z.extend(z_q11)

    print "Done."
    return z

'''
Creates fetal data for an entire missing chromosome from one parent
Final data includes all reads from one parent, and none from the other
'''
def longd(to_del, parent):
    print "Simulating complete deletion"
    z = filter(lambda x: x[3] != parent and x[3] != "-", to_del)
    print "Done."
    return z

'''
Returns fetal data for no aneuploidy
Final data includes all reads from both parents
'''
def noAneuploidy(fetal):
    print "Simulating no aneuploidy..."
    print "Done."
    return fetal

'''
Given plasma data, creates an observed sequence. For each position, sees if the number of reads is greater or less than the expected amount
Put into buckets of length 10,000 each
'''
def getSequence(data, ff):
    print "Generating observed sequence..."
    
    num_reads_p = len(data)*ff/2 # Paternal reads are 1/2 of the fetus'
    num_reads_m = len(data) - num_reads_p # Maternal reads are the rest
    expected_coverage_p = num_reads_p * READ_LEN * BUCKET_SIZE / GEN_LEN
    expected_coverage_m = num_reads_m * READ_LEN * BUCKET_SIZE / GEN_LEN
    
    # Generate two distributions, one for the father, and one for the mother
    low_p, high_p = poisson.interval(0.333, expected_coverage_p)
    low_m, high_m = poisson.interval(0.333, expected_coverage_m)
    
    # Count the number of times a read comes from m or p in each bucket
    coverage_p = defaultdict(0)
    coverage_m = defaultdict(0)
    for read in data:
        pos = int(read[0])
        read_len = len(read[1])
        for i in range(read_len):
            bucket = (pos + i)/BUCKET_SIZE
            if read[3] == "p":
                coverage_p[bucket] += 1
            else:
                coverage_m[bucket] += 1
        
    # Decide if the number of reads represents a low, normal, or high distribution     
    coverage = {}
    for key in coverage_p:
        if coverage_p[key] < low_p:
            p_val = "L"
        elif coverage_p[key] < high_p:
            p_val = "N"
        else:
            p_val = "H"
        if coverage_m[key] < low_m:
            m_val = "L"
        elif coverage_m[key] < high_m:
            m_val = "N"
        else:
            m_val = "H"
        val = (p_val, m_val)
        coverage[key] = val
        
    # Sort it by position
    observed_seq = []
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
    
def hammingDist(seq1, seq2):
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    # Make them the same distances
    if seq2_len > seq1_len:
        seq2 = seq2[:seq1_len]
    elif seq1_len > seq2_len:
        seq1 = seq1[:seq2_len]
    # Compute hamming distance
    dist = sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))
    return dist
    
def getDist(pos, seq, ref):
    pos = int(pos)
    low = 0
    high = len(ref) - 1
    while low <= high:
        mid = (low + high) / 2
        if int(ref[mid][0]) > pos:
            high = mid - 1
        elif int(ref[mid][0]) < pos:
            low = mid + 1
        else:
            # Match found - Compare the sequences
            dist = hammingDist(ref[mid][1], seq)
            break
    # A direct positional match has not been found
    if low > high:
        ref_seq = ref[mid][1]        
        diff = pos - int(ref[mid][0])
        if diff <= len(seq):
            seq_cpy = seq[:-diff]
            ref_cpy = ref_seq[diff:]
            dist = hammingDist(seq_cpy, ref_cpy)
        diff = int(ref[mid][0]) - pos
        if diff <= len(ref_seq):
            seq_cpy = seq[diff:]
            ref_cpy = ref_seq[:diff]
            dist2 = hammingDist(seq_cpy, ref_cpy)
    return min(dist, dist2)
    
def callReads(reads, m, p):
    print "Calling child reads..."
    for read in reads:
        pos = read[0]
        seq = read[1]
        dist_m = getDist(pos, seq, m)
        dist_p = getDist(pos, seq, p)
        if dist_m > dist_p:
            call = "p"
        elif dist_p > dist_m:
            call = "m"
        else:
            call = "-"
        read = read + (call,)
    print "Done."
    return reads

def main(ff, type, parent, display):
    m,f,p = loadGenomes()
    fetal = list(f)
    fetal = [tuple(l[0:-2].split(",")) for l in fetal]
    maternal = list(m)
    paternal = list(p)
    maternal = [tuple(l[0:-2].split(",")) + ("m",) for l in maternal]
    paternal = [tuple(l[0:-2].split(",")) + ("p",) for l in paternal]
    
    fetal = callReads(fetal, maternal, paternal)
    
    # Generate the aneuploidy for the entire fetus
    if type == "22q11del":
        fetal = del22q11(fetal, parent)
    elif type == "22q11dup":
        fetal = dup22q11(fetal, parent)
    elif type == "22q13del":
        fetal = del22q13(fetal, parent)
    elif type == "complete":
        fetal = complete(fetal, parent)
    elif type == "longd":
        fetal = longd(fetal, parent)
    elif type == "none":
        fetal = noAneuploidy(fetal)
    
    # Get the fetus samples that will show up in the plasma
    fetal = random.sample(fetal, int(READS*ff))
    
    # Fill the rest of the plasma reads
    g = copy.deepcopy(fetal)
    
    # Maternal DNA
    sample_size = READS - len(g)
    maternal_sample = random.sample(maternal, sample_size)
    g.extend(maternal_sample)
    
    #while(len(g) < READS):
    #    g.append(tuple(m.readline()[0:-2].split(",")) + ("m",))
    
    #sort the data on read position
    sorted(g, key=itemgetter(0))
    
    # Get an observation sequence
    seq = getSequence(g, ff)
    
    print "Writing to output file " + type + "_" + parent
    with open(OUTPUTPATH + type + "_" + parent, "w") as f:
        for entry in seq:
            f.write("%s,%s\n" % (entry[0], entry[1]))
    print "Done."
    
    if display:
        displayCoverage(g, type)


if __name__ == "__main__":
    args = parser.parse_args()
    main(normdist(FFMEAN,FFSTD)**2, args.t, args.p, args.d)
    #cProfile.run('main(normdist(FFMEAN,FFSTD)**2 * 100)')

