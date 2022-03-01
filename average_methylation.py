# input: list of samples, directory containing VCFs for SNPs for those samples and containing methylation results, window of how many bp surrounding those SNPs
#  output: VCF with added column for average percent methylation of methylation events in that window around those SNPS

import os
import gzip
import re
import sys

het_vcf_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/min_call_filter/01_process_dna/vqsr/"
meth_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/"
XX_placenta_samples = ["MW-11","MW-21","MW-31","OBG0055-P1"]
sampleid_to_patientid = {"MW-11":"OBG0088_Placenta", "MW-21":"OBG0095_Placenta","MW-31":"OBG0083_Placenta","OBG0055-P1":"OBG0055_Placenta"}
window = 5000 # look for methylation events within this many bp of the het SNP

#for sample in XX_placenta_samples:
sample = sys.argv[1]

# find het vcf in this directory for the selected samples
print ("Getting het vcf files in " + het_vcf_directory)
het_vcf_files = list()
chromosomes_found = list()
for (dirpath, dirnames, filenames) in os.walk(het_vcf_directory):
    het_vcf_files += [os.path.join(dirpath, file) for file in filenames if ".het.vcf" in file and not ".idx" in file and sample in file]

if (len(het_vcf_files) == 0):
    print ("Didn't find any het SNP vcf files here")
else: 
    print("Found: \n" + "\n".join(het_vcf_files))

for het_vcf_file in het_vcf_files: 
    print ("Processing " + het_vcf_file)

    # open output file
    outputfile = het_vcf_file.replace(".het.vcf",".het.meth.vcf")
    outputhandle = open(outputfile, "w")

    # read het SNP data
    het_SNP_handle = open (het_vcf_file, "r")
    het_SNP_data = het_SNP_handle.readlines()
    het_SNP_handle.close()

    # parse out the chromosome in file name

    find_chrom = re.search(r"chr\w", het_vcf_file)
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

    # for each het SNP, find the average % methylated for all methylation sites within the window specified

    position_index = 1  # HARD-CODED
    counter = 0
    for het_SNP in het_SNP_data:
        if (not "#" in het_SNP):
            SNP_items = het_SNP.replace("\n","").split("\t")
            pos = int(SNP_items[position_index])
            range_start = 0
            range_end = 0
            if (pos > window):
                range_start = pos - window
            range_end = pos + window

            meth_in_range_percentages = []
            percentage_index = 3  #HARD-CODED
            for meth_event in meth_data:
                meth_items = meth_event.replace("\n","").split("\t")
                if (int(meth_items[position_index]) >= range_start and int(meth_items[position_index]) <= range_end):
                     meth_in_range_percentages.append(float(meth_items[percentage_index]))
             
            if (len(meth_in_range_percentages) > 0):
                SNP_items.append(str(sum(meth_in_range_percentages)/len(meth_in_range_percentages)))
                print ("\t".join(SNP_items), file=outputhandle)
            else: 
                SNP_items.append("no_meth_info")
                print ("\t",join(SNP_items), file=outputhandle)

            counter += 1
            if (counter % 100 == 0):
                print ("Processed " + str(counter) + " out of " + str(len(het_SNP_data)))

    outputhandle.close()
    print ("Output created: " + outputfile)

