'''
Purpose: examine list of placenta specific genes in our RNAseq data,
        find the genes that have high expression in our placenta samples,
        look at their expression in decidua
Steps: 
1) for each gene in input list of published placenta specific genes
2) input directory of featurecounts results, grab out gene counts
3) look to see which of the genes in the placenta list are in the first column of the first gene counts file
4) for genes that are found, grab their data across all samples, with placenta samples first and decidua samples second
5) loop through genes and plot their values as a bar graph, different colors for placenta and decidua
6) for genes that are consistent in placenta, look to see which decidua have more expression, if decidua have consistent expression
7) estimate overall expression level in placenta by calculating range of average expression for all genes and calculating percentile of our gene, in first file
'''

import pandas as pd
import numpy as np 
import os

# load placenta specific gene list
placenta_gene_list_file = "/home/splaisie/scripts/placenta_gene_list.txt"
placenta_handle = open(placenta_gene_list_file, "r")
placenta_gene_list = placenta_handle.readlines()
placenta_handle.close()
placenta_gene_list = [x.replace("\n","") for x in placenta_gene_list]

# load the gene counts files from featurecounts results
placenta_featurecounts_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/feature_counts_rna/sccREF/"
placenta_fc_files = os.listdir(placenta_featurecounts_directory)
placenta_fc_genecounts_files = [placenta_featurecounts_directory + x for x in placenta_fc_files if ("geneCounts" in x and not "summary" in x and not "Decidua" in x)]

decidua_featurecounts_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/02_run_asereadcounter/feature_counts_rna/sccREF/"
decidua_fc_files = os.listdir(decidua_featurecounts_directory)
decidua_fc_genecounts_files = [decidua_featurecounts_directory + x for x in decidua_fc_files if ("geneCounts" in x and not "summary" in x and not ("Decidua" in x and "XY" in x))]

fc_genecounts_files = placenta_fc_genecounts_files + decidua_fc_genecounts_files

# in first placenta gene counts file in list
# a) look for which placenta genes are measured in our data

firstfile = fc_genecounts_files[0]
countdata = pd.read_table(firstfile, sep="\t",comment="#")
#countdata.rename({countdata.columns[6]:"Count"},axis=1,inplace=True)
countdata.drop(["Length","Chr","Start","End","Strand"],axis=1,inplace=True)

placenta_genes_found = []
for index,row in countdata.iterrows():
    if (row['Geneid'] in placenta_gene_list):
        placenta_genes_found.append(row['Geneid'])

# merge counts into one dataframe 
for fc_file in fc_genecounts_files[1:len(fc_genecounts_files)]:
    data_to_add = pd.read_table(fc_file, sep="\t",comment="#")
    data_to_add.drop(["Length","Chr","Start","End","Strand"],axis=1,inplace=True)
    countdata = pd.merge(countdata,data_to_add)

countdata.to_csv("/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/feature_counts_rna/sccREF/search_placenta_gene/countdata_merge.csv")

summary = countdata.describe()
summary.to_csv("/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/feature_counts_rna/sccREF/search_placenta_gene/countdata_summary.csv")

'''
for placenta_gene in placenta_genes_found:
    filt_countdata = countdata[countdata['Geneid'] == placenta_gene]

    ax = filt_countdata.plot(legend=False,kind="bar",logy=True,rot=90)
    fig = ax.get_figure()
    fig.savefig("/home/splaisie/scripts/placenta_genes_RNAseq." + placenta_gene + ".pdf")
'''
# load range of average gene expression in placenta samples
