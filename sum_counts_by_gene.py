# input: directory with _gene.tsv for annotated asereadcounter output
# ouput: _gene.tsv files with sums of ref counts, alt counts, and total counts, with a new allele balance of sum alt count / sum total counts
# steps: 
#   - for each _gene.tsv
#   - load into pandas data frame
#   - group by gene name
#   - take sums
#   - delete duplicates
#   - write to output file with _genesum.tsv ending

import os
import pandas as pd
import statistics as stats
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

inputdirectory = "/scratch/splaisie/placenta/valleywise/asereadcounter/analyze_ase_results/chrX/"
#inputfile = "/scratch/splaisie/placenta/valleywise/asereadcounter/analyze_ase_results/chrX/Plac_CON02_VW-ASUPlace-Cont02_180_013_S184_L001_chrX_allele_balance_gene_GRCh38.tsv"
genome = "GRCh38"

# get annotated tsv files
annotated_files = [inputdirectory+fn for fn in os.listdir(inputdirectory) if "gene_"+genome+".tsv" in fn]

for annotated_file in annotated_files: 
    outputfile = annotated_file.replace("_gene_"+genome+".tsv","_genesum_"+genome+".tsv")

    df = pd.read_csv(annotated_file, sep = "\t")
    df['sum_altcount'] = df.groupby(['Gene'])['alt_count'].transform('sum')
    df['sum_refcount'] = df.groupby(['Gene'])['ref_count'].transform('sum')
    df['sum_totalcount'] = df.groupby(['Gene'])['total_count'].transform('sum')

    new_df = df.drop_duplicates(subset=['Gene'])
    
    new_df['sum_max_altref'] = new_df[['sum_altcount','sum_refcount']].values.max(1)
    new_df['sum_allele_balance'] = new_df['sum_max_altref'] / new_df['sum_totalcount']

    new_df.to_csv(outputfile,sep="\t",index=False)

