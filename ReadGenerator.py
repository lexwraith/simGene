'''
Created on Mar 27, 2015
'''

from subprocess import call
from os.path import isfile
import os
from datetime import datetime

def removeFiles():
    if isfile("output_dat.aln"):
        os.remove("output_dat.aln")
    if isfile("output_dat.fq"):
        os.remove("output_dat.fq")
        
def extractReads():
    if not isfile("output_dat.sam"):
        return False
    filename = datetime.now()
    with open("output_dat.sam") as input_file:
        with open(filename, "a") as output_file:
            for line in input_file:
                #skip header files
                if line[0] == "@":
                    continue
                line_components = line.split()
                if len(line_components) < 4:
                    continue
                position = line_components[3]
                read = line_components[9]
                output_file.write(position + "    " + read)
    os.remove("output_dat.sam")

if __name__ == '__main__':
    i = 0
    while i < 100:
        return_code = call("art_illumina -sam -i hs_alt_CHM1_1.1_chr22.fa -l 25 -ss HS25 -f 10 -o output_dat", shell=True)
        if return_code != 0:
            print "Error occurred generating the reads"
            continue
        
        removeFiles()
        result = extractReads()
        if not result:
            print "Error occurred extracting the reads"
            continue
        i += 1