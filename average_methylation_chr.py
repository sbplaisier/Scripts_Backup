'''
Input: directory for methlyation results, sample 
Output: violin plot of percent methylation across events by chromosome
Steps:
1) for each sample
2) load methylation data
3) set dict of lists for percent methylation observed on each chromosome
4) calculate the average methylation of all events seen in each chromosome
6) when finished, make violin plot
    
'''

import os
import gzip
import re
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats

meth_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/"

sampleid_to_patientid = {"MW-11":"OBG0088_Placenta", "MW-21":"OBG0095_Placenta","MW-31":"OBG0083_Placenta","OBG0055-P1":"OBG0055_Placenta"}

#sample = "OBG0055-P1"
#sample = "MW-11"
#sample = "MW-21"
#sample = "MW-31"
sample_sex = "XX"
samples = sampleid_to_patientid.keys()

for sample in samples:
    print ("Processing " + sample)

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

    # read methylation data, divide into chromosomes
    include_chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]

    meth_file = meth_files[0] 
    meth_handle = gzip.open(meth_file,"rt")
    chr_meth_data = {} #key = chr__, value = list of average methylated
    chr_index = 0  #HARD-CODED
    percentage_index = 3  #HARD-CODED

    for line in meth_handle:
        meth_items = line.replace("\n","").split("\t")
        if (meth_items[chr_index] in include_chromosomes):
            if (meth_items[chr_index] in chr_meth_data):
                chr_meth_data[meth_items[chr_index]].append(float(meth_items[percentage_index]))
            else: 
                chr_meth_data[meth_items[chr_index]] = [float(meth_items[percentage_index])]
             
    meth_handle.close()
    print ("Loaded methylation file: " + meth_file)

    outputfile = meth_directory + "mannwhitneyu_autosome_chrX_" + sample + ".txt"
    outputhandle = open(outputfile,"w")
    for autosome in include_chromosomes:
        if (not autosome == "chrX"):
            print ("Comparing " + autosome)
            print (autosome, file = outputhandle)
            print(stats.mannwhitneyu(x=chr_meth_data[autosome],y=chr_meth_data["chrX"],alternative= 'two-sided'), file = outputhandle)
    outputhandle.close()

    # find the average % methylated for all methylation sites within each chromosome
    '''
    chr_meth = [] # list of chromosomes to x labels
    avg_meth = [] # avg percent methylation

    for c in chr_meth_data.keys():
        chr_meth.append(c)
        avg_perc = sum(chr_meth_data[c])/len(chr_meth_data[c])
        avg_meth.append(avg_perc)
    '''

    # make violin plot of average percent methylated in each chromosome

    '''
    plt.figure(figsize = (12,3))
    index_list = [*range(1, len(chr_meth_data.keys())+1, 1)] 
    fig, ax = plt.subplots()

    ax.set(title="Chromosomes: " + sample + "(" + sample_sex + ")")
    ax.violinplot(chr_meth_data.values(),
                showmeans=True,
                showmedians=False,
                showextrema=False,
                widths=0.75)
    ax.set_xticks(index_list)
    ax.set_xticklabels(chr_meth_data.keys(), rotation = 90)
    ax.set_ylabel("Average percent methylated")

    # report output 
    outputfile = meth_directory + "percent_methylation_chr_" + sample + ".png"
    plt.savefig(outputfile, bbox_inches='tight')

    print ("Output created: " + outputfile)
    '''
