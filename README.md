# simGene
A project for a Computational Biology class; Simulates and detects different types of aneuploidy on Chromosome 22 using an HMM and copy number variation.

General workflow is as follows:
1. Convert BAM to SAM files - one for each of the mother, father, and child
2. From each SAM file, filter to only include sequence, location, and chromosome location (i.e. q11)
3. Call the child's reads to either the father or mother
4. Simulate a specific type of aneuploidy using the child's called reads
5. Mix a fetal fraction of this data and the mother's to simulate maternal plasma
6. Generate an observed sequence, based on position, and the number of reads called to the mother/father, using buckets.
7. Train each HMM on its appropriate data
8. Run each test data on all HMMs, and accept that with highest probability

##Requirements
This project was developed on a high-capacity server with a 1TB SSD and 60GB of RAM. The SSD helped speed up computation, and while all 1TB was not needed, a large amount of disk space in general is. The .bam files are each ~9GB, and the associated .sam files are ~50GB. Enough space is needed to store three of these temporarily, before filtering. 
The large amount of RAM is more of a necessity for the project's successful running. Because we had the resources, we tended to just load the entire files into memory where possible. A machine with significantly less memory may not be able to run this successfully. 

##ReadFilter.py
Usage: python ReadFilter.py <input_file> <output_file> <filtered_output>
Given an input .sam file, this module filters it to only include read, sequence, and location information per read mapped to Chromosome 22. This module goes through input_file line by line, filtering anything not needed. There are two output files, output_file and filtered_output. output_file contains reads from all chromosomes, while filtered_output contains only reads for chromosome 22. output_file is kept only as a precaution. 

##ReadSampler.py
Usage: python ReadSampler.py <output> [-t <type>] [-p <parent>] [-d]
From filtered read files, it generates an observed sequence to represent the read coverage of each location in the chromosome (using buckets). This observed sequence is written to the output file specified by the output argument. The -t flag is used to indicate the type of aneuploidy to simulate, with a default of "none". The options supported are "none", "22q11del", "22q11dup", "2q13del", "complete", and "longd". The -p flag is used to specify which copy of the chromosome this aneuploidy is on. Valid options are "p" and "m", with the default being "p". The -d flag is used to indicate that the read coverage should be plotted and saved as an image. The file name will be saved as <output><type>.png. 
There are a number of important features in this module. First, it calls all child reads to either the mother or father, if it hasn't already. This is done by using a simple hamming distance on the sequence observed in the parents. If no sequence is available for that position, it uses a portion of the read available to get an estimate. 
It then simulates the specified aneuploidy type, given the type and parent arguments. For example, for 22q11del, it removes all reads from the parent specified that occur on q11. A fraction of this data is then combined with a fraction of the maternal data, to create a representation of maternal plasma.
The module then creates an observed sequence by counting the number of reads called to the mother/father in each bucket along the chromosome. These counts are then compared against a threshold (from the expected number of reads) to determine if they are one of three classifications - low, normal, or high. The combination of these for the mother and father in each bucket is made into one element of an observed sequence, which is the final output of the module.

##Analysis.py
Usage: python Analysis.py
This module reads all generated observed sequences, and uses them to train a set of 11 HMMs using Baum-Welch. Each HMM has two states - aneuploidy or not. The module then outputs in a csv file the probabilities of each observed sequence coming from each HMM in an aneuploidy state, using the forward algorithm. 
