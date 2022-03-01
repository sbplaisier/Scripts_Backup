# input: 
#   1) two directories containing .het.vcf files for heterozygous variants in individual samples
#   2) field you want to plot
#   3) chromosome to use (ex: chrX)
#  output: violin plot of that feature across variants in both batches of samples
#  steps:
#   1) read all lines of each file (minus those that start with #)
#   2) parse feature of interest
#   3) store list for each sample
#   4) output violin plot with all ranges found

import os
import matplotlib.pyplot as plt
import argparse

batch1directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/01_process_dna/vqsr/"
#batch2directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/vqsr/"
batch2directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/lesstrim/vqsr/"

#comparechr = "chrX"
parser = argparse.ArgumentParser(description = "violin plot comparing vcf feature between two batches of samples")
parser.add_argument('chr', help="Enter the chromosome to plot variants from (ex: chrX)")

#plotfeature = "BaseQRankSum"
parser.add_argument('feature', help="Enter the vcf feature to plot (ex: DP)")
args = parser.parse_args()

comparechr = args.chr
plotfeature = args.feature 

#outputfile = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/variant_list_compare_" + comparechr + "_" + plotfeature + ".pdf"
outputfile = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/lesstrim/variant_list_compare_lesstrim_" + comparechr + "_" + plotfeature + ".pdf"

def pull_feature_data(filedirectory,featurename, addlist,addlistnames):
    batchfilesall = os.listdir(filedirectory)
    batchfiles = [file for file in batchfilesall if "called.vqsr.sv.biallelic.snp" in file and ".het.vcf" in file and comparechr in file and not ".idx" in file]


    for file in batchfiles:
        variant_data = []
        inputhandle = open(filedirectory + file,"r")
        inputdata = inputhandle.readlines();
        inputhandle.close()

        for input in inputdata:
            if (not input[0] == "#"):
                input = input.replace("\n","")
                items = input.split("\t")
                variantinfo = items[7].split(";")
                for feature in variantinfo:
                    if (feature.find(featurename) == 0):
                        variant_data.append(float(feature.replace(featurename+"=","",1)))
        addlist.append(variant_data)
        short_name = file.replace(".gatk.called.vqsr.sv.biallelic.snp", "")
        short_name = short_name.replace(".het.vcf", "")
        addlistnames.append(short_name)


all_feature_data = list()
all_names = list()

pull_feature_data(batch1directory, plotfeature, all_feature_data, all_names)
count_batch1 = len(all_feature_data)
pull_feature_data(batch2directory, plotfeature, all_feature_data, all_names)

fig, ax = plt.subplots(figsize=(6, 3))

index_list = [*range(1, len(all_names)+1, 1)]
ax.set_xticklabels(all_names, rotation = 90)

violin_parts = ax.violinplot(all_feature_data,
            showmeans=False,
            showmedians=False,
            showextrema=False)
ax.set_xticks(index_list)
ax.set_xticklabels(all_names, rotation = 90)
ax.set_ylabel(plotfeature)

for pc in violin_parts['bodies']:
    pc.set_facecolor('darkblue')

plt.axvline(x = count_batch1+0.5, color = 'grey', linewidth=1)

plt.savefig(outputfile, bbox_inches='tight')
print("Plot output: " + outputfile)
