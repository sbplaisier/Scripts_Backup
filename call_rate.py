# input: 
#   1) sv.biallelic.vcf.gz for variants you have allele balance for
#   2) directory of allele balance info
#  output: tsv with another column added for call rate
#  steps:
#   1) read variant vcf and store number of times you find "./." in line by position
#   2) read through allele balance files
#   3) get position of each variant
#   4) print the call rate at the end

import os
import gzip

#vcf_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/vqsr/chr8.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
vcf_file = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/01_process_dna/vqsr/chr8.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
vcf_handle = gzip.open(vcf_file,"rt")
vcf_data = vcf_handle.readlines()
vcf_handle.close()

# dictionaries to hold the stuff about each variant to print out
call_rate = {}
qual = {}
info = {}
for vcf_line in vcf_data:
    if (not "#" in vcf_line):
        items = vcf_line.split ("\t")
        position = items[1]
        var_qual = items[5]
        var_info = items[7]

        no_call_count = vcf_line.count('./.')
        call_rate[position] = no_call_count

        qual[position] = var_qual
        info[position] = var_info

#ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/asereadcounter/analyze_ase_results/"
ase_directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/03_analyze_ase/analyze_ase_results/"
autosome = "chr8"

allele_balance_files = list()
for (dirpath, dirnames, filenames) in os.walk(ase_directory + autosome):
    allele_balance_files += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file and autosome in file]

if (len(allele_balance_files) == 0):
    print ("Didn't find any allele balance files here")

for allele_balance_file in allele_balance_files:
    inputhandle = open (allele_balance_file, "r")
    lines = inputhandle.readlines()
    inputhandle.close()
    
    outputfile = allele_balance_file.replace(".tsv","_callrate.tsv")
    outputhandle = open (outputfile, "w")

    print ("\t".join([lines[0].replace("\n",""), "Num_no_call","QUAL","INFO"]), file=outputhandle)
    
    for line in lines[1:len(lines)]:
        items = line.split("\t")
        position = items[1]
        if (not position in info):
            print ("Position " + position + " not found in vqsr variant file")
            continue
        else: 
            print ("\t".join([line.replace("\n",""), str(call_rate[position]), qual[position], info[position]]), file=outputhandle)

    outputhandle.close()
