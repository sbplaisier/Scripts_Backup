# input: 
#   - chrX or chr8
#   - nonSCC alignment bam filtered by chromosome entered
#   - SCC alignment bam filtered by chromosome entered
#   - entered these as parameters so we can script through the two reference genomes and chrX and chr8
# ouput: 
#   - table of reads with position and difference of alignment position of nonSCC alignment vs sex aware alignment
# steps: 
#   - make hash with key as the read sequence and value of position for nonSCC and SCC filtered alignment files
#   - for all reads in common, print the position in two alignments and calculated diff
#   - for reads that were only aligned in one or the alignment, set diff to a high number
#   - output table of values and scatter plot of position vs diff

import os
import pandas as pd
import statistics as stats
import numpy as np
import argparse

print ("Getting input")

chromosome_entered = "chrX"

#sex_aware_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.GRCh38.p12.genome.XX.sorted.chrX.txt"
#default_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.GRCh38.p12.genome.nonSCC.sorted.chrX.txt"

sex_aware_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.T2T.genome.XX.sorted.chrX.txt"
default_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.T2T.genome.nonSCC.sorted.chrX.txt"

#sex_aware_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.GRCh38.p12.genome.XX.sorted.chr8.txt"
#default_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.GRCh38.p12.genome.nonSCC.sorted.chr8.txt"

#sex_aware_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.T2T.genome.XX.sorted.chr8.txt"
#default_alignment = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/TCGA_LIHC_ED-A7PX-10A-01D-A34Z.T2T.genome.nonSCC.sorted.chr8.txt"

print ("Getting position of reads from alignments")

# alignment position hashes
sex_aware_position = {}
default_position = {}
print ("First, for SCC ref...")
sex_aware_input = open(sex_aware_alignment,"r")

for line in sex_aware_input:
    if (line[0] == '#'):
        continue

    line_split = line.split("\t")
    chr_found = line_split[2]
    position_found = line_split[3]
    readseq_found = line_split[9]

    #print (chr_found)
    #print(position_found)
    #print(readseq_found)

    if (not readseq_found in sex_aware_position):
        sex_aware_position[readseq_found] = position_found

sex_aware_input.close()

print ("Then, for default ref...")

default_input = open(default_alignment,"r")

for line in default_input:
    if (line[0] == '#'):
        continue

    line_split = line.split("\t")
    chr_found = line_split[2]
    position_found = line_split[3]
    readseq_found = line_split[9]

    #print (chr_found)
    #print(position_found)
    #print(readseq_found)
 
    if (not readseq_found in default_position):
        default_position[readseq_found] = position_found

default_input.close()

print ("Getting intersection of reads")

# intersection of hashes 

sex_aware_set = set(sex_aware_position)
default_set = set(default_position)

intersecting_reads = sex_aware_set.intersection(default_set)
sex_aware_specific = sex_aware_set.difference(default_set)
default_specific = default_set.difference(sex_aware_set)

#print (intersecting_reads.pop())
#print (sex_aware_specific.pop())
#print (default_specific.pop())

print ("Calculating differences between alignment positions for reads mapped in SCC and default alignments")

# base name for output file
outputfile = default_alignment.replace(".non_SCC","")
outputfile = outputfile.replace(".txt",".diffalign.csv")

df = pd.DataFrame(columns=['chromosome','read_sequence','SCC_position',"def_position","diff_def_minus_SCC"])

total = len(intersecting_reads)
count = 1
for read in intersecting_reads:
    SCC_pos_found = int(sex_aware_position[read])
    def_pos_found = int(default_position[read])
    diff_def_minus_SCC = def_pos_found - SCC_pos_found

    df = df.append({'chromosome':chromosome_entered,'read_sequence':read,'SCC_position':SCC_pos_found,"def_position":def_pos_found,"diff_def_minus_SCC":diff_def_minus_SCC}, ignore_index=True)
    if (count % 50000 == 0):
        print ("Processing " + str(count) + " of " + str(total))
        
        # print data frame in parts
        outputfile2 = outputfile.replace(".csv","_SCCspec"+ str(count) + ".csv")
        df.to_csv(outputfile2,sep="\t",index=False)
        df = pd.DataFrame(columns=['chromosome','read_sequence','SCC_position',"def_position","diff_def_minus_SCC"])

    count += 1

outputfile2 = outputfile.replace(".csv","_fileLast.csv")
df.to_csv(outputfile2,sep="\t",index=False)


print ("Reporting reads mapped in SCC alignment alone")

df = pd.DataFrame(columns=['chromosome','read_sequence','SCC_position',"def_position","diff_def_minus_SCC"])

total = len(sex_aware_specific)
count = 1
for read in sex_aware_specific:
    SCC_pos_found = int(sex_aware_position[read])
    def_pos_found = "not_found"
    diff_def_minus_SCC = -100000000 # randomchosen based on how big the differences in alignment were

    df = df.append({'chromosome':chromosome_entered,'read_sequence':read,'SCC_position':SCC_pos_found,"def_position":def_pos_found,"diff_def_minus_SCC":diff_def_minus_SCC}, ignore_index=True)
    if (count % 50000 == 0):
        print ("Processing " + str(count) + " of " + str(total))
        
        # print data frame in parts
        outputfile2 = outputfile.replace(".csv","_SCCspec"+ str(count) + ".csv")
        df.to_csv(outputfile2,sep="\t",index=False)
        df = pd.DataFrame(columns=['chromosome','read_sequence','SCC_position',"def_position","diff_def_minus_SCC"])

    count += 1

outputfile2 = outputfile.replace(".csv","_SCCspecLast.csv")
df.to_csv(outputfile2,sep="\t",index=False)

print ("Reporting reads mapped in default alignmnet alone")

df = pd.DataFrame(columns=['chromosome','read_sequence','SCC_position',"def_position","diff_def_minus_SCC"])

total = len(default_specific)
count = 1
for read in default_specific:
    SCC_pos_found = "not_found"
    def_pos_found = int(default_position[read])
    diff_def_minus_SCC = 100000000 # randomly chosen based on how big the differences in alignment were

    df = df.append({'chromosome':chromosome_entered,'read_sequence':read,'SCC_position':SCC_pos_found,"def_position":def_pos_found,"diff_def_minus_SCC":diff_def_minus_SCC}, ignore_index=True)
    if (count % 50000 == 0):
        print ("Processing " + str(count) + " of " + str(total))
        
        # print data frame in parts
        outputfile2 = outputfile.replace(".csv","_DEFspec"+ str(count) + ".csv")
        df.to_csv(outputfile2,sep="\t",index=False)
        df = pd.DataFrame(columns=['chromosome','read_sequence','SCC_position',"def_position","diff_def_minus_SCC"])

    count += 1

outputfile2 = outputfile.replace(".csv","_DEFspecLast.csv")
df.to_csv(outputfile2,sep="\t",index=False)


