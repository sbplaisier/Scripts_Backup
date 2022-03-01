# input: directory with unphased allele balance data (allele_balance.tsv from 03_analyze_ase in Placenta_XCI)
# output: scatter plot of allele balance across the chromosomes
# steps: 
# - read all allele balance files for chrX
# - find companion file for chr8
# - plot chromosome position on x-axis
# -    and "allele_balance" on y-axis
# - plot chrX and chr8 side by side for all samples

import os
import matplotlib.pyplot as plt

# input parameters

#ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/lesstrim/asereadcounter/analyze_ase_results/"
ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/min_call_filter/01_process_dna/asereadcounter/analyze_ase_results/"
autosome = "chr8"

#samples = ["MW-11","MW-21","MW-31","OBG0055-P1"]   #exome file ids
samples = ["MW-11","MW-21","MW-31","OBG0055-P1", "MW-53", "MW-43","MW-15","MW-33","OBG0055-D5"]   #exome file ids

# read allele balances 

allele_balance_files_chrX = list()
for exomesample in samples:
    for (dirpath, dirnames, filenames) in os.walk(ase_directory + "chrX/"):
        allele_balance_files_chrX += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file and exomesample in file and "chrX" in file]

for file in allele_balance_files_chrX:

    # parse out sample name
    fileitems = file.split("/")
    filenameitems = fileitems[-1].split("_")
    samplename = "_".join(filenameitems[0:3])

    # read allele balance from chrX file
    allele_balance_data_chrX = []
    position_data_chrX = []
    inputhandle = open(file,"r")
    inputdata = inputhandle.readlines();
    inputhandle.close()
    inputheaders = inputdata[0].replace('\n','').split('\t')
    allele_balance_index_chrX = inputheaders.index('allele_balance')
    position_index_chrX = inputheaders.index('position')

    for input in inputdata[1:len(inputdata)]:
        input = input.replace("\n","")
        items = input.split("\t")
        allele_balance_data_chrX.append(float(items[allele_balance_index_chrX]))
        position_data_chrX.append(int(items[position_index_chrX])/1000000)
   

    # read allele balance from chromosome input to compare to
    comparefile = file.replace("chrX",autosome)
    allele_balance_data_autosome = []
    position_data_autosome = []
    inputhandle = open(comparefile,"r")
    inputdata = inputhandle.readlines();
    inputhandle.close()
    inputheaders = inputdata[0].replace('\n','').split('\t')
    allele_balance_index_autosome = inputheaders.index('allele_balance')
    position_index_autosome = inputheaders.index('position')

    for input in inputdata[1:len(inputdata)]:
        input = input.replace("\n","")
        items = input.split("\t")
        allele_balance_data_autosome.append(float(items[allele_balance_index_autosome]))
        position_data_autosome.append(int(items[position_index_autosome])/1000000)
   

    # plot allele balance across samples

    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(12, 8))

    ax1.scatter(x = position_data_chrX, y = allele_balance_data_chrX, marker = 'o', s = 14, c = 'blue')
    
    ax1.set_ylabel("Unphased Allele balance")
    ax1.set_ylim(0,1)
    ax1.set_xlabel("chromosomeX position (Mb)")
    ax1.set_title ("Sample " + samplename)

    ax2.scatter(x = position_data_autosome, y = allele_balance_data_autosome, marker = 'o', s = 14, c = 'blue')

    ax2.set_ylabel("Unphased Allele balance")
    ax2.set_ylim(0,1)
    ax2.set_xlabel(autosome.replace("chr","chromosome") + " position (Mb)")
    ax2.set_title ("Sample " + samplename)

    # save figure
    
    fig.tight_layout(pad=2.0)
    outputpng = ase_directory + "unphased_allele_balance_acrosschr_" + samplename + ".png"
    plt.savefig(outputpng)
    print ("Created plot: ", outputpng)
