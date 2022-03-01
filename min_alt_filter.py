# input: 
#   1) vcf file you want filtered 
#   2) min thresold for filtering count of alt allele
#  output: 
#   1) filtered vcf file with min alt allele count at least threshold
#  steps:
#   1) for each line of vcf, read index 9 and on (genotypes)
#   2) keep count of how many 1's you see in genotypes
#   3) if >= threshold, print line, else skip

import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', type=str, required=True)
parser.add_argument('--threshold', type=int, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


#vcf_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/vqsr/chr8.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
#vcf_file = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/01_process_dna/vqsr/chr8.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"

# read lines of vcf file
vcf_data = []
vcf_file = args.vcf
if (".gz" in vcf_file):
    vcf_handle = gzip.open(vcf_file,"rt")
    vcf_data = vcf_handle.readlines()
    vcf_handle.close()
elif (".vcf" in vcf_file):
    vcf_handle =  open(vcf_file,"r")
    vcf_data = vcf_handle.readlines()
    vcf_handle.close()

# set outputfile
'''
if (".vcf.gz" in vcf_file):
    outputfile = vcf_file.replace(".vcf.gz",".minalt.vcf")
elif (".vcf" in vcf_file):
    outputfile = vcf_file.replace(".vcf",".minalt.vcf")
else: 
    print ("Input file not a vcf file")
    quit()
'''
outputfile = args.output
outputhandle = open (outputfile, "w")

# read vcf line and filter for min alt count below threshold
for vcf_line in vcf_data:
    if ("#" in vcf_line):
        print (vcf_line.replace("\n",""), file=outputhandle)
    elif (not "#" in vcf_line):
        items = vcf_line.split ("\t")
        
        alt_count = 0
        for genotype in items[9:len(items)]:
            gt_items = genotype.split(":")

            # count 1 in genotypes
            gt = gt_items[0]
            if "/" in gt:
                alleles = gt.split("/")
            elif "|" in gt:
                alleles = gt.split("|")
            for allele in alleles:
                if (allele == "1"):
                    alt_count += 1

        if (alt_count >= args.threshold):
            print (vcf_line.replace("\n",""), file=outputhandle)
        # else: print nothing and go to the next line

outputhandle.close()
