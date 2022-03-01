# input: chromosome you want to compare to chrX, starting directory with asereadcounter results,
#        directory with exon intersected bed files, list of samples
# output: scatter plots of allele balance in all samples across chr8 and chrX
# steps: 
# - for each sample
# - read and store chromosome positions of all the exon filtered variants 
#        for that sample
# - read all allele balance files in chrX subdirectory
# - store list of allele_balance
# - create list of point colors based on if the variant for which we 
#         measure allele specific expression is in the exon filtered
#         list or not
# - when all are collected, plot and save

import os
import matplotlib.pyplot as plt

# input parameters

autosome = "chr8"

#ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/asereadcounter/analyze_ase_results/"
ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/lesstrim/asereadcounter/analyze_ase_results/"
#ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/01_process_dna/asereadcounter/analyze_ase_results/"
#ase_directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/03_analyze_ase/analyze_ase_results/"
#outputpng = ase_directory + "allele_balance_exon_plot_" + autosome + ".png"
outputpng = ase_directory + "allele_balance_exon_plot_filtAltCount0_" + autosome + ".png"
#outputpng = ase_directory + "allele_balance_transcript_plot_" + autosome + ".png"

#vcf_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/vqsr/"
vcf_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/lesstrim/vqsr/"
#vcf_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/01_process_dna/vqsr/"
#vcf_directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/01_process_dna/vqsr/"

samples = ["MW-11","MW-21","MW-31","OBG0055-P1"]   #exome file ids
#samples = ["MW-15","MW-43","MW-33","OBG0055-D5","MW-53"]   #exome file ids
#samples = [
#"OBG0022",
#"OBG-0024-PLAC",
#"OBG-0026-PLAC",
#"OBG-0028-PLAC",
#"OBG-0030-PLAC",
#"OBG-0039-PLAC",
#"OBG0044",
#"OBG-0050-PLAC",
#"OBG-0051-PLAC",
#"OBG-0066-PLAC",
#"OBG0068",
#"OBG0111",
#"OBG0115",
#"OBG0120",
#"OBG-0121-PLAC",
#"OBG0133",
#"OBG-0138-PLAC",
#"OBG0156",
#"OBG0166",
#"OBG0170",
#"OBG0174",
#"OBG0175",
#"OBG0178",
#"OBG-0180-PLAC",
#"OBG-0188-PLAC",
#"OBG-0201-PLAC",
#"OBG-0205-PLAC",
#"OBG-0289-PLAC",
#"OBG-0338-PLAC",
#"OBG-0342-PLAC"
#]

# find the het.exon.bed files that contain the exon filtered variants for each sample
exonfiles_chrX = []
exonfiles_autosome = []
for exomesample in samples:
    for (dirpath, dirnames, filenames) in os.walk(vcf_directory):
        exonfiles_chrX += [os.path.join(dirpath, file) for file in filenames if "chrX" in file and "gatk.called.vqsr.sv.biallelic.snp" in file and "het.exon.bed" in file and exomesample in file]
#        exonfiles_chrX += [os.path.join(dirpath, file) for file in filenames if "chrX" in file and "gatk.called.vqsr.sv.biallelic.snp" in file and "het.transcript.bed" in file and exomesample in file]
        exonfiles_autosome += [os.path.join(dirpath, file) for file in filenames if autosome in file and "gatk.called.vqsr.sv.biallelic.snp" in file and "het.exon.bed" in file and exomesample in file]
#        exonfiles_autosome += [os.path.join(dirpath, file) for file in filenames if autosome in file and "gatk.called.vqsr.sv.biallelic.snp" in file and "het.transcript.bed" in file and exomesample in file]

# store the variants that are in exons (vqsr het.vcf bedtools interasected with annotation.exon)
exon_variants_chrX = {} # key = sample, value = list of position of chrX variants filtered by exons)
exon_variants_autosome = {} # key = sample, value = list of position of autosome variants filtered by exons)

# --- first for chrX
sample_index = 0 # since files were selected in order of the samples, keep track of index so you can look them up again
for exonfile_chrX in exonfiles_chrX:
    store_filtered_variants = []
    inputhandle = open(exonfile_chrX,"r")
    lines = inputhandle.readlines()
    inputhandle.close()
    for line in lines:
        items = line.split("\t")
        store_filtered_variants.append(items[2])
    exon_variants_chrX[samples[sample_index]] = store_filtered_variants
    sample_index += 1

# --- then for selected autosome
sample_index = 0 
for exonfile_autosome in exonfiles_autosome:
    store_filtered_variants = []
    inputhandle = open(exonfile_autosome,"r")
    lines = inputhandle.readlines()
    inputhandle.close()
    for line in lines:
        items = line.split("\t")
        store_filtered_variants.append(items[2])
    exon_variants_autosome[samples[sample_index]] = store_filtered_variants
    sample_index += 1


# read allele balances 

allele_balance_files_chrX = list()
for exomesample in samples:
    for (dirpath, dirnames, filenames) in os.walk(ase_directory + "chrX/"):
        allele_balance_files_chrX += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file and exomesample in file]

all_allele_balance = list()     # list of lists containing all allele balances given for each sample
all_names = list()              # list of sample names for which we collect allele balance
all_allele_balance_colors = list()        # list of lists containing point colors marking exons: dark red for exon, grey for not in exon

for file in allele_balance_files_chrX:

# --- parse out sample name
    fileitems = file.split("/")
    filenameitems = fileitems[-1].split("_")
    samplename = filenameitems[0]

# --- read allele balance from chrX file
    allele_balance_data = []
    allele_balance_colors = []
    inputhandle = open(file,"r")
    inputdata = inputhandle.readlines();
    inputhandle.close()
    inputheaders = inputdata[0].replace('\n','').split('\t')
    allele_balance_index = inputheaders.index('allele_balance')
    position_index = inputheaders.index('position')

    for input in inputdata[1:len(inputdata)]:
        input = input.replace("\n","")
        items = input.split("\t")
        allele_balance_data.append(float(items[allele_balance_index]))
        if (items[position_index] in exon_variants_chrX[samplename]):
            allele_balance_colors.append('blue')
        else:
            allele_balance_colors.append('grey')
   
    all_allele_balance.append(allele_balance_data)
    all_allele_balance_colors.append(allele_balance_colors)
    storename = file.replace(ase_directory + "chrX/","")
    storename = storename.replace("_allele_balance.tsv","")
    all_names.append(storename)

# --- read allele balance from chromosome input to compare to
    comparefile = file.replace("chrX",autosome)
    allele_balance_data = []
    allele_balance_colors = []
    inputhandle2 = open(comparefile,"r")
    inputdata2 = inputhandle2.readlines();
    inputhandle2.close()
    inputheaders = inputdata2[0].replace('\n','').split('\t')
    allele_balance_index = inputheaders.index('allele_balance')
    position_index = inputheaders.index('position')

    for input in inputdata2[1:len(inputdata2)]:
        input = input.replace("\n","")
        items = input.split("\t")
        allele_balance_data.append(float(items[allele_balance_index]))
        if (items[position_index] in exon_variants_autosome[samplename]):
            allele_balance_colors.append('blue') 
        else:
            allele_balance_colors.append('grey') 

    
    all_allele_balance.append(allele_balance_data)
    all_allele_balance_colors.append(allele_balance_colors)
    storename = comparefile.replace(ase_directory + autosome + "/","")
    storename = storename.replace("_allele_balance.tsv","")
    all_names.append(storename)


# plot allele balance across samples


fig, ax = plt.subplots(figsize=(16, 3))

ax.violinplot(all_allele_balance, 
            showmeans=False,
            showmedians=False,
            showextrema=False)

#ax.boxplot(all_allele_balance,
#            showmeans=False,
#            showfliers=False)

#create list of x and y coordinates to plot
xcoordinates = []
ycoordinates = []
pointcolors = []

for i in range(0,len(all_allele_balance)):
    for j in range(0,len(all_allele_balance[i])):
        #xcoordinates.append(i+1)
        ycoordinates.append(all_allele_balance[i][j])
        pointcolors.append(all_allele_balance_colors[i][j])
        if(all_allele_balance_colors[i][j] == 'blue'):
            xcoordinates.append(i+1.2)
        else:
            xcoordinates.append(i+1.1)

ax.scatter(x = xcoordinates, y = ycoordinates, marker = '.', s = 6, c = pointcolors)

index_list = [*range(1, len(all_names)+1, 1)]
ax.set_xticks(index_list)
ax.set_xticklabels(all_names, rotation = 90)
ax.set_ylabel("Unphased Allele balance")
ax.set_ylim(0,1)

plt.savefig(outputpng, bbox_inches='tight')
print ("Created plot: ", outputpng)
