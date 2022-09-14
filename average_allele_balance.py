# input: directory with unphased allele balance data (allele_balance.tsv from 03_anaalyze_ase in Placenta_XCI
# ouput: table of all the files with their mean, median, min, max, Q1, and Q3 for allele balances
# steps: 
# - get all allele_balance.tsv files
# - load table with pandas
# - get stats of allele balance column
# - print output

import os
import pandas as pd
import statistics as stats
import numpy as np

ase_directory = "/scratch/splaisie/placenta/valleywise/asereadcounter/analyze_ase_results/"
autosome = "chr8"

samples =  ["Plac_CON02", "Plac_CON03","Plac_CON05","Plac_CON06","Plac_CON10","Plac_HDP01","Plac_HDP08","Plac_HDP09","Plac_HDP10"]

# filter flag set to 1 if you want to filter by alt count
filter_alt_count = 1
alt_count_min_threshold = 2

# read allele balances 

allele_balance_files_chrX = list()
for exomesample in samples:
    for (dirpath, dirnames, filenames) in os.walk(ase_directory + "chrX/"):
        allele_balance_files_chrX += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file and exomesample in file and "chrX" in file]


# set output file
outputfile = ase_directory+"allele_balance_stats.tsv"
if (filter_alt_count):
    outputfile = outputfile.replace(".tsv","_filteraltcount.tsv")
outputhandle = open(outputfile, "w")
row_elements = ["file","chr","mean","median","min","max","Q1","Q3"]
print("\t".join(row_elements), file=outputhandle)

# iterate through chrX files and matching chr8 files

for chrXfile in allele_balance_files_chrX: 
    # get stats for allele balance for sites in chrX
    df = pd.read_csv(chrXfile, sep = "\t")
    data = df["allele_balance"]
    if (filter_alt_count):
        df2 = df.loc[df["alt_count"] >= alt_count_min_threshold]
        data = df2["allele_balance"]
    
    meanAB = stats.mean(data)
    medianAB = stats.median(data)
    minAB = min(data)
    maxAB = max(data)
    Q1_AB = np.percentile(data,25)
    Q3_AB = np.percentile(data,75)
    row_elements = [chrXfile.replace(ase_directory,""),"chrX",str(meanAB),str(medianAB),str(minAB),str(maxAB),str(Q1_AB),str(Q3_AB)]
    print("\t".join(row_elements),file=outputhandle)

    # get stats for allele balance for sites in chr8
    chr8file = chrXfile.replace("chrX","chr8")
    df = pd.read_csv(chr8file, sep = "\t")
    data = df["allele_balance"]
    if (filter_alt_count):
        df2 = df.loc[df["alt_count"] >= alt_count_min_threshold]
        data = df2["allele_balance"]

    meanAB = stats.mean(data)
    medianAB = stats.median(data)
    minAB = min(data)
    maxAB = max(data)
    Q1_AB = np.percentile(data,25)
    Q3_AB = np.percentile(data,75)
    row_elements = [chr8file.replace(ase_directory,""),"chr8",str(meanAB),str(medianAB),str(minAB),str(maxAB),str(Q1_AB),str(Q3_AB)]
    print("\t".join(row_elements),file=outputhandle)

