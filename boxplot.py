import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import ttest_ind

inputdir = "C:\\Users\\splaisie\\Dropbox (ASU)\\PC\\Documents\\Placenta\\asereadcounter\\analyze_ase_results\\"
chrXfiles = os.listdir(inputdir + "chrX\\")

for file in chrXfiles: 
    #inputfile = "chrX\\MW-15_OBG0088_Decidua2_chrX_allele_balance.tsv"
    #inputfile2 = "chr8\\MW-15_OBG0088_Decidua2_chr8_allele_balance.tsv"
    #title = "MW-15_OBG0088_Decidua2"
    
    chrXfile = inputdir + "chrX\\" + file
    chr8file = inputdir + "chr8\\" + file.replace("chrX","chr8")
    title = file.replace("_chrX_allele_balance.tsv","")
    
    df = pd.read_csv(chrXfile, sep="\t")
    df2 = pd.read_csv(chr8file, sep="\t")
    chrX_allele_balance_data = df['allele_balance']
    chr8_allele_balance_data = df2['allele_balance']
    t, p = ttest_ind(chrX_allele_balance_data, chr8_allele_balance_data)
    
    plotdata = []
    plotdata.append(chrX_allele_balance_data)
    plotdata.append(chr8_allele_balance_data)
    
    fig = plt.figure(figsize = (4,6))
    
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xticklabels(["chrX","chr8"])
    ax.set_ylabel("allele_balance")
    ax.set_title(title + " (ttest p = " + "{:.2e}".format(p) + ")")
    
    bp = ax.boxplot(plotdata)
    
    plt.savefig(inputdir + title + ".png", bbox_inches='tight')