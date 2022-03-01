# input: 
#   1) directory of allele balance info
#  output: 
#   1) file with gene feature that position is in 
#       plus other info to annotate that position
#  steps:
#   1) load the annotation file
#   2) for each line, get the chromosome and position
#   3) find entries that have that position
#   4) report 1, prioritizing exons, UTR, intron, then other

import os
import re

#ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/asereadcounter/phased_allele_balance/"
#ase_directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/04_phasing/phased_allele_balance/"
ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/min_call_filter/01_process_dna/asereadcounter/phased_allele_balance/"

# find allele balance files in this directory
print ("Getting allele balance files")
allele_balance_files = list()
chromosomes_found = list()
for (dirpath, dirnames, filenames) in os.walk(ase_directory):
    allele_balance_files += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file or "allele_balance_data.tsv" in file]

if (len(allele_balance_files) == 0):
    print ("Didn't find any allele balance files here")


# find which chromosomes we have allele balance files for

for file in allele_balance_files:    
    find_chrom = re.search(r"chr\w", file)
    if (not find_chrom.group() in chromosomes_found):
        chromosomes_found.append(find_chrom.group())
print ("Found data for chromosomes: " + ",".join(chromosomes_found))

# check to see what the delimiter of the files is
# sometimes file is labeled tsv when it is really a csv
inputhandle = open (allele_balance_files[0], "r")
lines = inputhandle.readlines()
inputhandle.close()
headerline = lines[0]
headers_tab = headerline.split("\t")
headers_comma = headerline.split(",")
delimiter = "\t"
if (len(headers_comma) > len(headers_tab)):
    delimiter = ","

# load annotation file

print ("Loading annotation file")
#annotation_file = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf"
annotation_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/annotated_genomes/gencode.v39.annotation.gtf"
annotation_handle = open(annotation_file,"r")
annotation_data = annotation_handle.readlines()
annotation_handle.close()

for current_chromosome in chromosomes_found:

    print ("Working on " + current_chromosome)
    # slice out annotations for the chromosomes we have data for
    filt_annotation = [line for line in annotation_data if line.split("\t")[0] == current_chromosome]

    # go through allele balance files annotating the sites

    for allele_balance_file in allele_balance_files:
        inputhandle = open (allele_balance_file, "r")
        lines = inputhandle.readlines()
        inputhandle.close()
        
        outputfile = allele_balance_file.replace(".tsv","_annotated_v39.tsv")
        outputhandle = open (outputfile, "w")

        headerline = lines[0].replace("\n","")
        headers = headerline.split(delimiter)
        headers.append("Feature")
        headers.append("Exon_number")
        headers.append("Gene_info")
        print ("\t".join(headers), file=outputhandle)
        
        for line in lines[1:len(lines)]:
            lineitems = line.replace("\n","").split(delimiter)
            position = lineitems[1]
            feature = lineitems[2]
            found = []  # set to feature types found
            chosen_info = ""
            for annotation_line in filt_annotation:
                items = annotation_line.split("\t")
                start = items[3]
                end = items[4]

                if (int(position) >= int(start) and int(position) <= int(end)):
                    if (not "CDS" in found or not "exon" in found):  # keep looking until you find a CDS or exon
                        feature = items[2]
                        info = items[8]
                        exon_num = ""
                        chosen_info = feature + "\t" 
                        if ("exon_number" in info):
                            #index_exon_num = info.find("exon_number ") + 12
                            #exon_num = info[index_exon_num] + "\t"
                            exon_search = re.search(r'exon_number \d+',info)
                            exon_num = exon_search.group() 
                        else: 
                            exon_num = "not found"
                        chosen_info += exon_num + "\t" + info.replace("\n","")  
                        found.append(feature)
                    
            if (found == False):
                lineitems.append("not found")
                lineitems.append("\t")
                lineitems.append("\t")
            else:
                lineitems.append(chosen_info)

            print ("\t".join(lineitems), file=outputhandle)

            
        outputhandle.close()

