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

A number of installation requirements are needed. The project runs assuming Python2.7. Also needed to be installed are:

1.  Samtools

2.  Numpy

3.  SciPy

4.  Matplotlib

5.  HMM from Colorado State (available here: http://www.cs.colostate.edu/~anderson/cs440/index.html/doku.php?id=notes:hmm2)

##Data
Original data is gotten from the publicly available 1000 genomes project. Three datasets are used: HG00689 (Father), HG00690 (Mother), and HG00691 (Child). These are available for download from the following paths (note: each file is ~9GB). 

HG00689: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00689/cg_data/HG00689_lcl_SRR822306.mapped.COMPLETE_GENOMICS.CGworkflow2_2_evidenceSupport.CHS.high_coverage.20130501.bam

HG00690: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00690/cg_data/HG00690_lcl_SRR822311.mapped.COMPLETE_GENOMICS.CGworkflow2_2_evidenceSupport.CHS.high_coverage.20130401.bam

HG00691: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00691/cg_data/HG00691_lcl_SRR826849.mapped.COMPLETE_GENOMICS.CGworkflow2_2_evidenceSupport.CHS.high_coverage.20130401.bam

##Full Working Example
Because of the very specific nature of our project, and the fact that it relies on dealing with big data, it does not make sense to run it on small sample inputs of data, so they have not been included.

What we have included is three files - child_filtered, father_filtered, and mother_filtered. These are a very small sample of the output of ReadFilter.py, and may be used to run ReadSampler. Additionally, we have included two sample output files. 22q11dupm0 is a sample observed sequence, generated for trisomy on 22q11 on the maternal side. We have also included heatmap.csv, the output of the entire project. 

Now, we will describe how to run the program on full data. In this example, we will assume that the downloaded .bam fiels are HG00689.bam, HG00690.bam, and HG00691.bam. 

samtools view -h -o father.sam HG00689.bam

python ReadFilter.py father.sam father father_filtered

samtools view -h -o mother.sam HG00690.bam

python ReadFilter.py mother.sam mother mother_filtered

samtools view -h -o child.sam HG00691.bam

python ReadFilter.py child.sam child child_filtered

python ReadSampler.py 22q11del_output -t 22q11del -p m

python Analysis.py


##ReadFilter.py
Usage: python ReadFilter.py <input_file> <output_file> <filtered_output>

Given an input .sam file, this module filters it to only include read, sequence, and location information per read mapped to Chromosome 22. This module goes through input_file line by line, filtering anything not needed. There are two output files, output_file and filtered_output. output_file contains reads from all chromosomes, while filtered_output contains only reads for chromosome 22. output_file is kept only as a precaution. 

##ReadSampler.py
Usage: python ReadSampler.py <output> [-t <type>] [-p <parent>] [-d]

From filtered read files, it generates an observed sequence to represent the read coverage of each location in the chromosome (using buckets). This observed sequence is written to the output file specified by the output argument. The -t flag is used to indicate the type of aneuploidy to simulate, with a default of "none". The options supported are "none", "22q11del", "22q11dup", "2q13del", "complete", and "longd". The -p flag is used to specify which copy of the chromosome this aneuploidy is on. Valid options are "p" and "m", with the default being "p". The -d flag is used to indicate that the read coverage should be plotted and saved as an image. The file name will be saved as <output>.png. 

There are a number of important features in this module. First, it calls all child reads to either the mother or father, if it hasn't already. This is done by using a simple hamming distance on the sequence observed in the parents. If no sequence is available for that position, it uses a portion of the read available to get an estimate. 

It then simulates the specified aneuploidy type, given the type and parent arguments. For example, for 22q11del, it removes all reads from the parent specified that occur on q11. A fraction of this data is then combined with a fraction of the maternal data, to create a representation of maternal plasma.

The module then creates an observed sequence by counting the number of reads called to the mother/father in each bucket along the chromosome. These counts are then compared against a threshold (from the expected number of reads) to determine if they are one of three classifications - low, normal, or high. The combination of these for the mother and father in each bucket is made into one element of an observed sequence, which is the final output of the module.

In order for correct connection with Analysis.py the output filename should reflect the type of the aneuploidy generated. This is of the form <type><parent><integer>, i.e. "22q11dupm5". 

##Analysis.py
Usage: python Analysis.py

This module reads all generated observed sequences in the output directory, and uses them to train a set of 11 HMMs using Baum-Welch. THe model to train is determined using filenames, so it is important to adhere to file naming conventions (i.e. 22q11dupm5). Each HMM has two states - aneuploidy or not. The module then outputs in a csv file the probabilities of each observed sequence coming from each HMM in an aneuploidy state, using the forward algorithm. The output of the program is heatmap.csv, written in the directory Analysis.py is run from.

##Config.py
This module contains all important magic numbers and file paths. For example, it defines where to programs should look for reads and write output (OUTPUTPATH). 

##Credit/hmm.py
The HMM module used was provided by Prof. Chuck Anderson of Colorado State University. It was received from his publicly-available site listed above. All credit belongs to him. 
