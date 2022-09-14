# input: directory with annotated and summed allele balance data (_gene_genome.tsv)
# ouput: list of genes with proportion of samples that had summed allele balance > 0.8
# steps: 
# - get all gene_genome.tsv files
# - get list of unique genes listed at least once in the gene_genome.tsv files 
# - load table with pandas
# - for each gene, iterate through each sample searching for the gene
# - in total number of normotensive, hypertensive, and all
# - count how many samples are allele balance > 0.8, median allele balance > 0.8
# - calculate proportions
# - write inactivated for >= 7/10 or 3/5 allele balance > 0.8, escape for <= 3/10 or 2/5, variable otherwise 
# - print output

import os
import pandas as pd
import statistics as stats
import numpy as np

ase_directory = "/scratch/splaisie/placenta/valleywise/asereadcounter/analyze_ase_results/"
autosome = "chr8"
xci_threshold = 0.8

samples =  ["Plac_CON02", "Plac_CON03","Plac_CON05","Plac_CON06","Plac_CON10","Plac_HDP01","Plac_HDP08","Plac_HDP09","Plac_HDP10"]

annotation = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf"
genome = "GRCh38"

chromosome_entered = "chrX"
gene_list_file = "/scratch/splaisie/placenta/valleywise/asereadcounter/analyze_ase_results/chrX/all_genes_removedups.txt"
gene_input = open(gene_list_file,"r")
gene_list = gene_input.readlines()
gene_list = [gene.replace("\n","") for gene in gene_list]

annotated_snp_files = [fn for fn in os.listdir(ase_directory+"chrX/")
              if any(sample in fn and "_genesum_"+genome+".tsv" in fn for sample in samples)]

out_df = pd.DataFrame(columns=['chromosome','gene','status_all',"status_norm","status_hyp",'medianAB_all',"medianAB_norm","medianAB_hyp"])

for gene in gene_list: 
    norm_AB = []
    hyp_AB = []
    all_AB = []
    for annotated_snp_file in annotated_snp_files:

        # read data for sites in chrX
        df = pd.read_csv(ase_directory+"chrX/"+annotated_snp_file, sep = "\t")
        df = df.loc[df["Gene"] == gene]
        if (len(df) > 0):
            all_AB.append(float(df["sum_allele_balance"]))

            if ("CON" in annotated_snp_file and "Cont" in annotated_snp_file):
                norm_AB.append(float(df["sum_allele_balance"]))
            elif ("HDP" in annotated_snp_file):
                hyp_AB.append(float(df["sum_allele_balance"]))
                

    median_norm_AB = 0.5
    status_norm_AB = "variable"
    if (len(norm_AB) > 0):
        median_norm_AB = stats.median(norm_AB)
        perchigh_norm_AB = len([element for element in norm_AB if element > xci_threshold]) / len(norm_AB)
        if (perchigh_norm_AB >= 0.7 and median_norm_AB >= xci_threshold):
            status_norm_AB = "inactivated"
        elif (perchigh_norm_AB <= 0.3 and median_norm_AB <= 0.75):
            status_norm_AB = "escape"

    median_hyp_AB = 0.5
    status_hyp_AB = "variable"
    if (len(hyp_AB) > 0):
        median_hyp_AB = stats.median(hyp_AB)
        perchigh_hyp_AB = len([element for element in hyp_AB if element > xci_threshold]) / len(hyp_AB)
        if (perchigh_hyp_AB >= 0.7 and median_hyp_AB >= xci_threshold):
            status_hyp_AB = "inactivated"
        elif (perchigh_hyp_AB <= 0.3 and median_hyp_AB <= 0.75):
            status_hyp_AB = "escape"

    median_all_AB = 0.5
    status_all_AB = "variable"
    if (len(all_AB) > 0):
        median_all_AB = stats.median(all_AB)
        perchigh_all_AB = len([element for element in all_AB if element > xci_threshold]) / len(all_AB)
        if (perchigh_all_AB >= 0.7 and median_all_AB >= xci_threshold):
            status_all_AB = "inactivated"
        elif (perchigh_all_AB <= 0.3 and median_all_AB <= 0.75):
            status_all_AB = "escape"

    out_df = out_df.append({'chromosome':chromosome_entered,'gene':gene,'status_all':status_all_AB,"status_norm":status_norm_AB,"status_hyp":status_hyp_AB,'medianAB_all':median_all_AB,"medianAB_norm":median_norm_AB,"medianAB_hyp":median_hyp_AB}, ignore_index=True)

    #print(all_AB)
    #print (out_df)


outputfile = gene_list_file.replace(".txt","_xci.txt")
out_df.to_csv(outputfile,sep="\t",index=False)
