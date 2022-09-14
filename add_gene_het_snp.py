# input: directory with unphased allele balance data (allele_balance.tsv from 03_anaalyze_ase in Placenta_XCI
# ouput: tsv file with the gene added in
# steps: 
# - get all allele_balance.tsv files
# - load table with pandas
# - go through all rows looking up what gene the SNP is in based on its position
# - print output

import os
import pandas as pd
import statistics as stats
import numpy as np

ase_directory = "/scratch/splaisie/placenta/valleywise/asereadcounter/analyze_ase_results/"
autosome = "chr8"

samples =  ["Plac_CON02", "Plac_CON03","Plac_CON05","Plac_CON06","Plac_CON10","Plac_HDP01","Plac_HDP08","Plac_HDP09","Plac_HDP10"]

annotation = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gtf"
genome = "GRCh38"

gene_lookup_chrX = {}
exon_lookup_chrX = {}
gene_lookup_autosome = {}
exon_lookup_autosome = {}

print ("Getting annotation  info")
with open(annotation,"r") as annotation_input:
    for line in annotation_input:
        if (not line.startswith("#")): 
            line_split = line.split("\t")
            chr_found = line_split[0]
            type_found = line_split[2]
            start_found = line_split[3]
            end_found = line_split[4]
            desc_found = line_split[8]
            desc_split = desc_found.split(";")
            
            if (type_found == "gene"):
                gene_name_found = desc_split[2]
                gene_name_found = gene_name_found.replace("gene_name ","")
                gene_name_found = gene_name_found.replace("\"","")
                gene_name_found = gene_name_found.replace(" ","")
            elif (type_found == "exon"):
                exon_name_found = desc_split[3]
                exon_name_found = exon_name_found.replace("gene_name ","")
                exon_name_found = exon_name_found.replace("\"","")
                exon_name_found = exon_name_found.replace(" ","")
                exon_number_found = desc_split[6]
                exon_number_found = exon_number_found.replace("exon_number ","")
                exon_number_found = exon_number_found.replace("\"","")
                exon_number_found = exon_number_found.replace(" ","")

            if(chr_found == "chrX"):
                if (type_found == "gene"):
                    if (not start_found in gene_lookup_chrX):
                        gene_lookup_chrX[start_found] = [end_found,gene_name_found]
                elif (type_found == "exon"):
                    if (not start_found in exon_lookup_chrX):
                        exon_lookup_chrX[start_found] = [end_found,exon_name_found,exon_number_found]
            if(chr_found == autosome):
                if (type_found == "gene"):
                    if (not start_found in gene_lookup_autosome):
                        gene_lookup_autosome[start_found] = [end_found,gene_name_found]
                elif (type_found == "exon"):
                    if (not start_found in exon_lookup_autosome):
                        exon_lookup_autosome[start_found] = [end_found,exon_name_found,exon_number_found]
                        

# filter flag set to 1 if you want to filter by alt count
filter_alt_count = False
alt_count_min_threshold = 2

chrXfiles = [fn for fn in os.listdir(ase_directory+"chrX/")
              if any(sample in fn and "allele_balance.tsv" in fn for sample in samples)]


for chrXfile in chrXfiles:

    # read data for sites in chrX

    print ("Annotating chrX")
    df = pd.read_csv(ase_directory+"chrX/"+chrXfile, sep = "\t")

    # filter on minimum alt_count value if parameters are set to do that
    if (filter_alt_count):
        df = df.loc[df["alt_count"] >= alt_count_min_threshold]

    # add blank column for gene name
    df['Gene'] = "not_found"
    df['GeneCoordinates'] = "not_found"
    df['ExonCoordinates'] = "not_found"

    # iterate through the 'position' column and look up the gene based on coordinates
    row_index = 0
    for position in df['position']:
        gene_matched = ""
        coordinates = ""
        exons = ""
        for start_pos in gene_lookup_chrX:
            if (int(position) >= int(start_pos) and int(position) <= int(gene_lookup_chrX[start_pos][0])):
                gene_matched = gene_lookup_chrX[start_pos][1]
                coordinates = start_pos + "-" + gene_lookup_chrX[start_pos][0]
                exon_gene_matched = exon_lookup_chrX[start_pos][1]
                exon_number_matched = exon_lookup_chrX[start_pos][2]
                exons +=  "(" + exon_gene_matched + " exon " +  exon_number_matched + ")"

        df.loc[row_index,"Gene"] = gene_matched
        df.loc[row_index,"GeneCoordinates"] = coordinates
        df.loc[row_index,"ExonCoordinates"] = exons
        row_index += 1

    outputfile = chrXfile.replace(".tsv","_gene_"+genome+".tsv")
    df.to_csv(ase_directory+"chrX/"+outputfile,sep="\t",index=False)
    
    # ---- do the same thing for autosome ----

    print ("Annotating chr8")
    # read data for sites in chr8
    chr8file = chrXfile.replace("chrX","chr8")
    df2 = pd.read_csv(ase_directory+"chr8/"+chr8file, sep = "\t")

    # filter on minimum alt_count value if parameters are set to do that
    if (filter_alt_count):
        df2 = df2.loc[df2["alt_count"] >= alt_count_min_threshold]

    # add blank column for gene name
    df2['Gene'] = "not_found"
    df2['GeneCoordinates'] = "not_found"
    df2['ExonCoordinates'] = "not_found"

    # iterate through the 'position' column and look up the gene based on coordinates
    row_index = 0
    for position in df2['position']:
        gene_matched = ""
        coordinates = ""
        exons = ""
        for start_pos in gene_lookup_autosome:
            if (int(position) >= int(start_pos) and int(position) <= int(gene_lookup_autosome[start_pos][0])):
                gene_matched = gene_lookup_autosome[start_pos][1]
                coordinates = start_pos + "-" + gene_lookup_autosome[start_pos][0]
                exon_gene_matched = exon_lookup_autosome[start_pos][1]
                exon_number_matched = exon_lookup_autosome[start_pos][2]
                exons +=  "(" + exon_gene_matched + " exon " +  exon_number_matched + ")"
        df2.loc[row_index,"Gene"] = gene_matched
        df2.loc[row_index,"GeneCoordinates"] = coordinates
        df.loc[row_index,"ExonCoordinates"] = exons
        row_index += 1

    outputfile = chr8file.replace(".tsv","_gene_"+genome+".tsv")
    df2.to_csv(ase_directory+"chr8/"+outputfile,sep="\t",index=False)
    

