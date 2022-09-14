'''
Input: directories containing XCI data files 
Output: violin plot showing average percent methylation across gene and promoter coordinates in samples
        broken up by XX and XY in placenta and decidua (on same range so we can compare)
Steps:
1) for each escape, inactivated, and variable (3 panel figure)
2) go through the directories loading data for the samples and selected genes
3) make violin plots for XX placenta, then XY placenta, decidua with XX placenta, with XY placenta
    
'''

import os
import gzip
import re
import sys
import matplotlib.pyplot as plt

selected_gene_file = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/Placenta_Xinactivation_Genes.csv"
meth_directories = ["/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/", "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/"]

sampleid_to_patientid = {"MW-11":"OBG0088_Placenta", "MW-21":"OBG0095_Placenta","MW-31":"OBG0083_Placenta","OBG0055-P1":"OBG0055_Placenta","MW-41":"OBG0023_Placenta","MW-55":"OBG0085_Placenta","MW-33":"OBG0083_Decidua", "MW-15":"OBG0088_Decidua","OBG0055-D5":"OBG0055_Decidua","MW-43":"OBG0023_Decidua","MW-53":"OBG0085_Decidua"}
placenta_XX_samples = ["MW-11", "MW-21","MW-31","OBG0055-P1"]
placenta_XY_samples = ["MW-41","MW-55"]
decidua_w_XX_samples = ["MW-33", "MW-15","OBG0055-D5"]
decidua_w_XY_samples = ["MW-43","MW-53"]

samples = sampleid_to_patientid.keys()

# go through directories finding XCI_data
xci_data_files = list()
for current_directory in meth_directories:
    for (dirpath, dirnames, filenames) in os.walk(current_directory):
        xci_data_files += [os.path.join(dirpath, file) for file in filenames if "XCI_data" in file and not "all" in file and not "wilcoxen" in file]

ordered_xci_data_files = list()
for sample in sampleid_to_patientid.keys():
    str_match = [s for s in xci_data_files if sample in s]
    ordered_xci_data_files.append(str_match[0])

gene_types = ["escape","inactivated","variable"]
for gt in gene_types: 
    gene_data_to_plot = []
    prom_data_to_plot = []
    sample_names = []

    print ("Processing " + gt)

    # get data for this gene type from the files
    for data_file in ordered_xci_data_files:
        gene_data = []
        promoter_data = []
        inputhandle = open(data_file,"r")
        lines = inputhandle.readlines()
        for line in lines:
            line = line.replace("\n","")
            if ("gene" in line and "XCI_status" in line and "avg_meth" in line):
                continue
            (gene_name, gene_status, gene_avg, prom_avg) = line.split("\t")
            if (gene_status == gt):
                gene_data.append(float(gene_avg))
                promoter_data.append(float(prom_avg))

        gene_data_to_plot.append(gene_data)
        prom_data_to_plot.append(promoter_data)
        file_split = data_file.split("/")
        file_name = file_split[-1]
        sample_name = file_name.replace("XCI_data_","")
        sample_name = sample_name.replace(".txt","")
        sample_names.append(sampleid_to_patientid[sample_name])


    # violin plot
    plt.figure(figsize = (6,3))
    fig, ax = plt.subplots(nrows=2, ncols=1)

    index_list = [*range(1, len(gene_data_to_plot)+1, 1)] 

    ax[0].set(title= gt + " genes")
    ax[0].violinplot(gene_data_to_plot,
                showmeans=False,
                showmedians=True,
                showextrema=False)
    #            showfliers=False)
    ax[0].set_xticks([])
    #ax[0].set_xticklabels(sample_names, rotation = 90)
    ax[0].set_ylabel("Average percent methylated")
    ax[0].set_ylim([0,100])
    ax[0].text(1.5,5,"XX placenta",fontsize = 8)
    ax[0].axvline(4.5)
    ax[0].text(4.75,5,"XY placenta",fontsize = 8)
    ax[0].axvline(6.5)
    ax[0].text(6.75,5,"Decidua with XX", fontsize = 8)
    ax[0].axvline(9.5)
    ax[0].text(9.6,5,"Decidua with XY", fontsize = 8)

    ax[1].set(title= gt + " promoters")
    ax[1].violinplot(prom_data_to_plot,
                showmeans=False,
                showmedians=True,
                showextrema=False)
    #            showfliers=False)
    ax[1].set_xticks(index_list)
    ax[1].set_xticklabels(sample_names, rotation = 90)
    ax[1].set_ylabel("Average percent methylated")
    ax[1].set_ylim([0,100])
    ax[1].text(1.5,5,"XX placenta",fontsize = 8)
    ax[1].axvline(4.5)
    ax[1].text(4.75,5,"XY placenta",fontsize = 8)
    ax[1].axvline(6.5)
    ax[1].text(6.75,5,"Decidua with XX", fontsize = 8)
    ax[1].axvline(9.5)
    ax[1].text(9.6,5,"Decidua with XY", fontsize = 8)

    # report output 
    outputfile = meth_directories[0] + "compare_meth_" + gt + "_genes.png"
    #outputfile = meth_directories[0] + "compare_meth_" + gt + "_genes_box.png"
    plt.savefig(outputfile, bbox_inches='tight')

    print ("Output created: " + outputfile)
