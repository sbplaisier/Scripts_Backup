# input: 
#   1) sv.biallelic.vcf.gz for variants you have allele balance for
#   2) directory of allele balance info
#  output: tsv with another column added for:
#   1) counts of number not called
#   2) num called 0 / (2 alleles x num samples)
#   3) num called 1 / (2 alleles x num samples)
#  steps:
#   1) read variant vcf and store number output fields by position
#   2) read through allele balance files
#   3) get position of each variant
#   4) print the additional info at the end

import os
import gzip

vcf_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/vqsr/chr8.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
#vcf_file = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/01_process_dna/vqsr/chr8.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
vcf_handle = gzip.open(vcf_file,"rt")
vcf_data = vcf_handle.readlines()
vcf_handle.close()

# dictionaries to hold the stuff about each variant to print out
qual = {}
info = {}
ref_allele_count = {} # count how many alleles in total genotypes are 0
alt_allele_count = {} # count how many alleles in total genotypes are 1
no_call_count = {} # count how many alleles in total genotypes are 1
AD_siteA_avg = {} # average AD (unfiltered allelic depth) first site across all genotypes
AD_siteB_avg = {} # average AD second site across all genotypes
DP_avg = {} # average DP (read depth at this position) across all genotypes
GQ_avg = {} # average GQ (conditional genotype quality) across all genotypes
for vcf_line in vcf_data:
    if (not "#" in vcf_line):
        items = vcf_line.split ("\t")
        position = items[1]
        var_qual = items[5]
        var_info = items[7]
        qual[position] = var_qual
        info[position] = var_info

        nocall_count = 0
        ref_count = 0
        alt_count = 0
        AD_siteA = []
        AD_siteB = []
        DP = []
        GQ = []
        for genotype in items[9:len(items)]:
            gt_items = genotype.split(":")

            # count 0, 1, . in genotypes
            gt = gt_items[0]
            if "/" in gt:
                alleles = gt.split("/")
            elif "|" in gt:
                alleles = gt.split("|")
            for allele in alleles:
                if (allele == "0"):
                    ref_count += 1
                elif (allele == "1"):
                    alt_count += 1
                elif (allele == '.'):
                    nocall_count += 1

            # calculate averages for other fields
            ad = gt_items[1]
            (adA,adB) = ad.split(',')
            AD_siteA.append(int(adA))
            AD_siteB.append(int(adB))

            dp = gt_items[2]
            if(not dp == '.'):
                DP.append(int(dp))

            gq = gt_items[3]
            if(not gq == '.'):
                GQ.append(int(gq))

        ref_allele_count[position] = ref_count
        alt_allele_count[position] = alt_count
        no_call_count[position] = nocall_count
        AD_siteA_avg[position] = sum(AD_siteA)/len(AD_siteA)
        AD_siteB_avg[position] = sum(AD_siteB)/len(AD_siteB)
        DP_avg[position] = sum(DP)/len(DP)
        GQ_avg[position] = sum(GQ)/len(GQ)


ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/asereadcounter/analyze_ase_results/"
#ase_directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/03_analyze_ase/analyze_ase_results/"
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
    
    outputfile = allele_balance_file.replace(".tsv","_allele_frequency.tsv")
    outputhandle = open (outputfile, "w")

    print ("\t".join([lines[0].replace("\n",""), "QUAL","INFO","Num_ref_call","Num_alt_call","Num_no_call", "avgAD_siteA","avgAD_siteB","avg_DP","avg_GQ"]), file=outputhandle)
    
    for line in lines[1:len(lines)]:
        items = line.split("\t")
        position = items[1]
        if (not position in info):
            print ("Position " + position + " not found in vqsr variant file")
            continue
        else: 
            print ("\t".join([line.replace("\n",""), qual[position], info[position], str(ref_allele_count[position]), str(alt_allele_count[position]),str(no_call_count[position]),str(AD_siteA_avg[position]), str(AD_siteB_avg[position]), str(DP_avg[position]), str(GQ_avg[position])]), file=outputhandle)

    outputhandle.close()
    quit()
