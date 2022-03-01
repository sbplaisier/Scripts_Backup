# input: file of genomic ranges, directory for methlyation results, samples to run
#  output: genomic range file with average percent methylation in events within each range tacked on

import os
import gzip
import re
import sys

range_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/min_call_filter/01_process_dna/asereadcounter/phased_allele_balance_filterPARs/promoters/OBG0055-P1_chrX_phased_allele_balance_data_annotated.gene_name.promoters.csv"
meth_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/"

sampleid_to_patientid = {"MW-11":"OBG0088_Placenta", "MW-21":"OBG0095_Placenta","MW-31":"OBG0083_Placenta","OBG0055-P1":"OBG0055_Placenta"}
window = 5000 # look for methylation events within this many bp of the het SNP

#for sample in XX_placenta_samples:
#sample = sys.argv[1]
sample = "OBG0055-P1"

print ("Processing " + range_file)

# open output file
outputfile = range_file.replace(".csv",".meth.csv")
outputhandle = open(outputfile, "w")

# read range data
range_handle = open (range_file, "r")
range_data = range_handle.readlines()
range_handle.close()

# parse out the chromosome in file name

find_chrom = re.search(r"chr\w", range_file)
chrom_found = find_chrom.group()
print ("Chromosome: " + chrom_found)

# load methylation file (bismark coverage .cov.gz)
meth_files = list()
for (dirpath, dirnames, filenames) in os.walk(meth_directory):
    meth_files += [os.path.join(dirpath, file) for file in filenames if "deduplicated.bismark.cov.gz" in file and sampleid_to_patientid[sample] in file]

if (len(meth_files) > 1):
    print ("Multiple methylation files found, not sure which to use: " + "\t".join(meth_files))
    quit()
elif (len(meth_files) == 0):
    print ("Didn't find any methylation (bismark coverage) files for the sample " + sample)
    quit()
else:
    print ("Matching methylation file: \n" + "\n".join(meth_files))

# read methylation data
meth_file = meth_files[0] 
meth_handle = gzip.open(meth_file,"rt")
#meth_data = meth_handle.readlines()
meth_data = []
for line in meth_handle:
    if (chrom_found in line):
        meth_data.append(line.replace("\n",""))
meth_handle.close()

print ("Loaded methylation file: " + meth_file)

# for each range, find the average % methylated for all methylation sites within that range

start_position_index = 1  # HARD-CODED
end_position_index = 2  # HARD-CODED
counter = 0
for range_line in range_data:
    if (not "#" in range_line):
        range_items = range_line.replace("\n","").split(",")
        range_start = int(range_items[start_position_index])
        range_end = int(range_items[end_position_index])

        meth_in_range_percentages = []
        position_index = 1 # HARD-CODED
        percentage_index = 3  #HARD-CODED
        for meth_event in meth_data:
            meth_items = meth_event.replace("\n","").split("\t")
            if (int(meth_items[position_index]) >= range_start and int(meth_items[position_index]) <= range_end):
                 meth_in_range_percentages.append(float(meth_items[percentage_index]))
         
        if (len(meth_in_range_percentages) > 0):
            range_items.append(str(sum(meth_in_range_percentages)/len(meth_in_range_percentages)))
            print (",".join(range_items), file=outputhandle)
        else: 
            range_items.append("no_meth_info")
            print (",",join(SNP_items), file=outputhandle)

        counter += 1
        if (counter % 100 == 0):
            print ("Processed " + str(counter) + " out of " + str(len(range_data)))

outputhandle.close()
print ("Output created: " + outputfile)

