import pandas as pd
import matplotlib.pyplot as plt
import os

inputdir = "C:\\Users\\splaisie\\Dropbox (ASU)\\PC\\Documents\\Placenta\\asereadcounter\\analyze_ase_results\\"
chrXfiles = os.listdir(inputdir + "chrX\\")

chrX_data = []
chr8_data = []
label_data = []
for file in chrXfiles: 
    if ("Placenta" in file):
        chrXfile = inputdir + "chrX\\" + file
        chr8file = inputdir + "chr8\\" + file.replace("chrX","chr8")
        title = file.replace("_chrX_allele_balance.tsv","")
        
        chrXmatrix = pd.read_csv(chrXfile, sep="\t")
        chrX_allele_balance_data = chrXmatrix['allele_balance']
        chrX_data.append(chrX_allele_balance_data)
        
        chr8matrix = pd.read_csv(chr8file, sep="\t")
        chr8_allele_balance_data = chr8matrix['allele_balance']
        chr8_data.append(chr8_allele_balance_data)
        
        label_data.append(title)
    
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(25, 4))

axs[0].violinplot(chrX_data,
                  showmeans=False,
                  showmedians=True,
                  showextrema=False)
axs[0].set_title('X Chromosome')


axs[1].violinplot(chr8_data,
                  showmeans=False,
                  showmedians=True,
                  showextrema=False)
axs[1].set_title('Chromosome 8')

xtick_locations = [y + 1 for y in range(len(label_data))]
for ax in axs:
    ax.yaxis.grid(True)
    ax.set_xticks(xtick_locations)
    ax.set_xticklabels(label_data, rotation = 90)
    ax.set_ylim (0,1.2)

plt.savefig(inputdir + "AlleleBalance_violin.png", bbox_inches='tight')