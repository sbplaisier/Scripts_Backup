# input: directory with alignment difference output (*diffalign*)
# ouput: image of scatter plot
# steps: 
# - load all the 'def_position' (SCC_position if not_found) and 'diff_def_minus_SCC' from the *diff_align*csv files into 2 long lists
# - scatter plot with matplotlib

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import collections as mc

data_directory = "/data/CEM/wilsonlab/projects/TCGA_LIHC/alignments/dna/"
sample = "ED-A7PX-10A-01D-A34Z"

diff_align_files = [fn for fn in os.listdir(data_directory) if sample in fn and "diffalign" in fn and "csv" in fn]

x = []
y = []

for diff_align_file in diff_align_files:
    df = pd.read_csv(data_directory+diff_align_file, sep = "\t")
    pos_data = df['def_position']
    if (pos_data[2] == "not_found"):
        pos_data = df['SCC_position']
    x = np.append(x,pos_data)

    diff_data = df['diff_def_minus_SCC']
    y = np.append(y,diff_data)

point_size = [3] * len(x)

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,gridspec_kw= {'height_ratios':[4,1]}, figsize=(8,5))

# add scatter plot of alignment differences
#ax1 = plt.subplot(211)
ax1.scatter(x, y, s=point_size)
ax1.set_title ("Alignment position differences on GRCh38")
ax1.set_xlabel("Read alignment position")
ax1.set_ylabel("difference (default - sex aware)")
#ax1.ticklabel_format(useOffset=False,style='plain')

# add regions of chrX

#ax2 = plt.subplot(212,sharex = ax1)
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/
chrX_region_lines = [[(10001, 1.2), (2781479, 1.2)], [(319338, 1.3), (601516, 1.3)], [(4950957, 1.3), (5129468, 1.3)], [(58605580, 1), (62412542, 1)], [(79965154, 1.3), (80097082, 1.3)], [(155701383, 1.2), (156030895, 1.2)]]
chrX_region_names = ["PAR_1","REGION_187","REGION_239","CENX","REGION_188","PAR_2"]
chrX_region_colors = np.array([(1, 0, 0, 1),(0, 1, 0, 1),(0, 1, 0, 1),(0, 0, 1, 1),(0, 1, 0, 1),(1, 0, 0, 1)])

lc = mc.LineCollection(chrX_region_lines, colors=chrX_region_colors,linewidths=15)
ax2.add_collection(lc)
ax2.set_ylim([0.9,1.4])
ax2.set_xticks([])
ax2.set_yticks([])

plt.savefig(data_directory+"diffalign_scatterplot_"+sample+".png")
