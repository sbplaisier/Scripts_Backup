# input: chromosome you want to compare to chrX, starting directory with asereadcounter results
#        intersected with transcript coordinates
# output: density plots (like violin plots) of allele balance in all samples
# steps: 
# 1) read all transcripts.rmdups.bed files in chrX subdirectory
# 2) iterate through, for each one
# 3) store list of allele_balance
# 4) find corresponding chr__ transcript.bed file
# 5) store list of allele balance there
# 6) when all are collected, plot and save

import os
import matplotlib.pyplot as plt

directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/lesstrim2/asereadcounter/analyze_ase_results/"
#directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/03_analyze_ase/analyze_ase_results/"
comparechr = "chr8"

samples = ["MW-11","MW-21","MW-31","OBG0055-P1"]   #exome file ids
#samples = ["OBG0044","OBG0068","OBG0111","OBG0115","OBG0120","OBG0133","OBG0156","OBG0166","OBG0170","OBG0174","OBG0175","OBG0178","OBG0022","OBG-0024-PLAC","OBG-0028-PLAC","OBG-0030-PLAC","OBG-0039-PLAC","OBG-0050-PLAC","OBG-0051-PLAC","OBG-0066-PLAC","OBG-0121-PLAC","OBG-0138-PLAC","OBG-0180-PLAC","OBG-0188-PLAC","OBG-0201-PLAC","OBG-0205-PLAC","OBG-0289-PLAC","OBG-0338-PLAC","OBG-0342-PLAC"]
 
outputpng = directory + "allele_balance_plot_" + comparechr + ".png"

chrXfiles = list()
for exomesample in samples:
    for (dirpath, dirnames, filenames) in os.walk(directory + "chrX/"):
        chrXfiles += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file and exomesample in file]

all_allele_balance = list()
all_names = list()
for file in chrXfiles:
    #print ("chrxfile " + file)

    # read allele balance from chrX file
    allele_balance_data = []
    inputhandle = open(file,"r")
    inputdata = inputhandle.readlines();
    inputhandle.close()

    for input in inputdata[1:len(inputdata)]:
        input = input.replace("\n","")
        items = input.split("\t")
        allele_balance_data.append(float(items[-1])) #assumption: allele balance is last column
   
    all_allele_balance.append(allele_balance_data)
    storename = file.replace(directory + "chrX/","")
    storename = storename.replace("_allele_balance.tsv","")
    all_names.append(storename)

    # read allele balance from chromosome input to compare to
    comparefile = file.replace("chrX",comparechr)
    #print (comparechr + "file " + comparefile)
    allele_balance_data = []
    inputhandle = open(comparefile,"r")
    inputdata = inputhandle.readlines();
    inputhandle.close()

    for input in inputdata[1:len(inputdata)]:
        input = input.replace("\n","")
        items = input.split("\t")
        allele_balance_data.append(float(items[-1])) #assumption: allele balance is last column
    
    all_allele_balance.append(allele_balance_data)
    storename = comparefile.replace(directory + comparechr + "/","")
    storename = storename.replace("_allele_balance.tsv","")
    all_names.append(storename)

#print(all_allele_balance)
#plt.boxplot(all_allele_balance,
#        showfliers=False)
#plt.xticks(all_names)

index_list = [*range(1, len(all_names)+1, 1)]

#fig, (ax1,ax2) = plt.subplots(2,1,figsize=(6, 3))
fig, ax = plt.subplots(figsize=(12, 3))

#ax1.boxplot(all_allele_balance,
#                  showmeans=False,
#                  showfliers=False)
#ax1.tick_params(axis="x",labelbottom="False")
#ax1.set_xticklabels([])
#ax1.set_ylabel("Allele balance")
#ax1.set_ylim(0,1)

ax.violinplot(all_allele_balance,
            showmeans=False,
            showmedians=False,
            showextrema=False)
ax.set_xticks(index_list)
ax.set_xticklabels(all_names, rotation = 90,fontsize=5)
ax.set_ylabel("Unphased Allele balance")
#a2.set_ylim(0,1)

plt.savefig(outputpng, bbox_inches='tight')

print("Created plot: " + outputpng)
