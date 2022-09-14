'''
Input: directory for methlyation results, samples to compare, chrX genes to run, annotation files for 
        genes and transcripts, num upstream bases for promoter
Output: violin plot showing average percent methylation across gene and promoter coordinates in samples
Steps:
1) load the coordinates of genes and promoters for chrX genes entered
    a) key = gene_name in chrX, value = tuple of (start,end), first entry in chrX
    b) gene hash: first entry in annotation.gene.gtf
    c) promoter hash: first entry in annotation.transcript.gtf, (start of transcript - num_upstream, start of transcript)
2) for each sample
3) load methylation data
4) for each gene
5) calculate the average methylation in the gene and its promotor
6) when finished, make panel of violin plots
    a) two violin plots for each sample 
    b) one by gene coordinates
    c) one by promoter coordinates
    
'''

import os
import gzip
import re
import sys
import matplotlib.pyplot as plt

selected_gene_file = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/Placenta_Xinactivation_Genes.csv"
meth_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/"
#meth_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/"
gene_annotation_file = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gene.gtf"
transcript_annotation_file = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.transcript.gtf"
num_upstream = 2000 # num bases upstream of transcription start site to represent promoter

sampleid_to_patientid = {"MW-11":"OBG0088_Placenta", "MW-21":"OBG0095_Placenta","MW-31":"OBG0083_Placenta","OBG0055-P1":"OBG0055_Placenta","MW-41":"OBG0023_Placenta","MW-55":"OBG0085_Placenta"}
XX_samples = ["MW-11", "MW-21","MW-31","OBG0055-P1"]
XY_samples = ["MW-41","MW-55"]
'''
sampleid_to_patientid = {"MW-33":"OBG0083_Decidua", "MW-15":"OBG0088_Decidua","OBG0055-D5":"OBG0055_Decidua","MW-43":"OBG0023_Decidua","MW-53":"OBG0085_Decidua"}
XX_samples = ["MW-33", "MW-15","OBG0055-D5"]
XY_samples = ["MW-43","MW-53"]
'''

#sample = "OBG0055-P1"
#sample = "MW-11"
#sample = "MW-21"
#sample = "MW-31"
#sample_sex = "XX"
samples = sampleid_to_patientid.keys()

# flags to indicate what output you want
outputDataFile = False
plotData = True 
#for sample in samples: 
for sample in XX_samples: 
    print ("Processing " + sample)
    if (sample in XX_samples):
        sample_sex = "XX"
    elif (sample in XY_samples):
        sample_sex = "XY"

    # load selected genes: chrX genes that are assigned to be Escape, Inactivated, or Variable_escape
    select_handle = open(selected_gene_file,"r")
    select_data = select_handle.readlines()
    select_handle.close()

    selected_genes = {}
    escape_genes = []
    xinact_genes = []
    varescape_genes = []
    for xin in select_data:
        xin = xin.replace("\n","")
        (name,status) = xin.split(",")
        if (status == "Escape"):
            escape_genes.append(name)
        elif (status == "Inactivated"):
            xinact_genes.append(name)
        elif (status == "Variable_escape"):
            varescape_genes.append(name)

    selected_genes = {}
    selected_genes["escape"] = escape_genes
    selected_genes["inactivated"] = xinact_genes
    selected_genes["variable"] = varescape_genes


    # collect gene list and indices needed to read and filter annotation gtf files
    all_selected_genes = selected_genes["escape"] + selected_genes["inactivated"] + selected_genes["variable"]

    chr_index = 0 # HARD-CODED for gtf files
    start_index = 3
    end_index = 4
    info_index = 8 
    chrom_selected = "chrX"

    # load gene coordinates for selected genes
    gene_handle = open(gene_annotation_file,"r")
    gene_data = gene_handle.readlines()
    gene_handle.close()

    gene_coordinates = {} # key = gene name, value = (start,end)
    for annotated_gene in gene_data:
        # get needed data from annotated gene
        annotated_gene = annotated_gene.replace("\n","")
        gene_items = annotated_gene.split("\t")

        chr_found = gene_items[chr_index]
        start_found = gene_items[start_index]
        end_found = gene_items[end_index]
        info_found = gene_items[info_index]

        # parse out gene name from gene info
        find_gene_name = re.search(r"gene_name \"\w+", info_found)
        gene_name_found = find_gene_name.group()
        gene_name_found = gene_name_found.replace("gene_name \"","")
        
        if (gene_name_found in all_selected_genes and chr_found == chrom_selected):
            gene_coordinates[gene_name_found] = (int(start_found), int(end_found))

    # load promotor coordinates for selected genes' transcripts
    transcript_handle = open(transcript_annotation_file,"r")
    transcript_data = transcript_handle.readlines()
    transcript_handle.close()

    promoter_coordinates = {} # key = gene name, value = (start - num_unstream,start)
    for annotated_transcript in transcript_data:
        # get needed data from annotated gene
        annotated_transcript = annotated_transcript.replace("\n","")
        gene_items = annotated_transcript.split("\t")

        chr_found = gene_items[chr_index]
        start_found = gene_items[start_index]
        end_found = gene_items[end_index]
        info_found = gene_items[info_index]

        # parse out gene name from gene info
        find_gene_name = re.search(r"gene_name \"\w+", info_found)
        gene_name_found = find_gene_name.group()
        gene_name_found = gene_name_found.replace("gene_name \"","")
        
        if (gene_name_found in all_selected_genes and chr_found == chrom_selected):
            promoter_coordinates[gene_name_found] = (int(start_found) - num_upstream, int(start_found))


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

    # read methylation data for selected chromosome
    meth_file = meth_files[0] 
    meth_handle = gzip.open(meth_file,"rt")
    meth_data = []
    for line in meth_handle:
        if (chrom_selected in line):
            meth_data.append(line.replace("\n",""))
    meth_handle.close()

    print ("Loaded methylation file: " + meth_file)


    # for each gene included in analysis
    # find the average % methylated for all methylation sites within that gene and its promoter

    avg_meth = {} # key = gene_name, value (avg percent methylation gene, avg perc methylation promotor)

    for g in all_selected_genes:
        if (g in gene_coordinates and g in promoter_coordinates):
            (gene_start,gene_end) = gene_coordinates[g]
            (promoter_start,promoter_end) = promoter_coordinates[g]

            gene_percentages = []
            promoter_percentages = []

            position_index = 1 # HARD-CODED
            percentage_index = 3  #HARD-CODED

            for meth_event in meth_data:
                meth_items = meth_event.replace("\n","").split("\t")
                if (int(meth_items[position_index]) >= gene_start and int(meth_items[position_index]) <= gene_end):
                     gene_percentages.append(float(meth_items[percentage_index]))
                if (int(meth_items[position_index]) >= promoter_start and int(meth_items[position_index]) <= promoter_end):
                     promoter_percentages.append(float(meth_items[percentage_index]))
             

            if (len(gene_percentages) > 0):
                avg_gene_perc = sum(gene_percentages)/len(gene_percentages)
            else: 
                avg_gene_perc = 0

            if (len(promoter_percentages) > 0):
                avg_promoter_perc = sum(promoter_percentages)/len(promoter_percentages)
            else: 
                avg_promoter_perc = 0

            avg_meth[g] = (avg_gene_perc, avg_promoter_perc)

    # for genes and promoters, get distributions of escape, inactivated, and variable_escape for violin plot and print to file

    if (outputDataFile):
        dataoutfile = meth_directory + "XCI_data_" + sample + ".txt"
        datahandle = open(dataoutfile,"w")
        print ("\t".join(["gene","XCI_status","gene_avg_meth","promoter_avg_meth"]), file= datahandle)

    escape_genes_meth = []
    inactive_genes_meth = []
    variable_genes_meth = []
    escape_promoters_meth = []
    inactive_promoters_meth = []
    variable_promoters_meth = []

    for g in avg_meth:
        if (g in selected_genes["escape"] and g in avg_meth):
            (gene_avg_meth,promoter_avg_meth) = avg_meth[g]
            escape_genes_meth.append(gene_avg_meth)
            escape_promoters_meth.append(promoter_avg_meth)
            if (outputDataFile):
                print ("\t".join([g,"escape",str(gene_avg_meth),str(promoter_avg_meth)]), file= datahandle)
        elif (g in selected_genes["inactivated"] and g in avg_meth):
            (gene_avg_meth,promoter_avg_meth) = avg_meth[g]
            inactive_genes_meth.append(gene_avg_meth)
            inactive_promoters_meth.append(promoter_avg_meth)
            if (outputDataFile):
                print ("\t".join([g,"inactivated",str(gene_avg_meth),str(promoter_avg_meth)]), file= datahandle)
        elif (g in selected_genes["variable"] and g in avg_meth):
            (gene_avg_meth,promoter_avg_meth) = avg_meth[g]
            variable_genes_meth.append(gene_avg_meth)
            variable_promoters_meth.append(promoter_avg_meth)
            if (outputDataFile):
                print ("\t".join([g,"variable",str(gene_avg_meth),str(promoter_avg_meth)]), file= datahandle)
            
    if (outputDataFile):
        datahandle.close()

    # make violin plot of average methylation in the different gene groups using gene coordinates and promoter coordinates

    if (plotData):
        gene_meth_data = list() 
        gene_meth_data.append(escape_genes_meth)
        gene_meth_data.append(inactive_genes_meth)
        gene_meth_data.append(variable_genes_meth)

        plt.figure(figsize = (6,3))
        index_list = [*range(1, len(selected_genes)+1, 1)] #len(selected_genes) = num categories of selected genes = 3 = escape, inactivated, variable
        #fig, ax = plt.subplots(nrows=1, ncols=2, sharey = 'row')
        fig, ax = plt.subplots(nrows=1, ncols=2)

        #ax[0].set(title="Genes:" + sampleid_to_patientid[sample] + "(" + sample_sex + ")")
        ax[0].set(title=sampleid_to_patientid[sample] + "(" + sample_sex + ")")
        ax[0].violinplot(gene_meth_data,
                    showmeans=True,
                    showmedians=False,
                    showextrema=False)
        ax[0].set_xticks(index_list)
        #ax[0].set_xticklabels(["Escape (n = 21)","Inactive (n = 132)","Variable_escape (n = 37)"], rotation = 90)
        ax[0].set_xticklabels(["ESC","XCI","VAR"], rotation = 90)
        ax[0].set_ylabel("Average percent methylated")

        promoter_meth_data = list() 
        promoter_meth_data.append(escape_promoters_meth)
        promoter_meth_data.append(inactive_promoters_meth)
        promoter_meth_data.append(variable_promoters_meth)

        #ax[1].set(title="Promotors:" + sampleid_to_patientid[sample] + "(" + sample_sex + ")")
        ax[1].set(title= sampleid_to_patientid[sample] + "(" + sample_sex + ")")
        ax[1].violinplot(promoter_meth_data,
                    showmeans=True,
                    showmedians=False,
                    showextrema=False)
        ax[1].set_xticks(index_list)
        #ax[1].set_xticklabels(["Escape (n = 21)","Inactive (n = 132)","Variable_escape (n = 37)"], rotation = 90)
        ax[1].set_xticklabels(["ESC","XCI","VAR"], rotation = 90)
        ax[1].set_ylabel("Average percent methylated")

        # report output 
        outputfile = meth_directory + "xinactivation_genes_showmeans_" + sample + ".png"
        plt.savefig(outputfile, bbox_inches='tight')

        print ("Output created: " + outputfile)

