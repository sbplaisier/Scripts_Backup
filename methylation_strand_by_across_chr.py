# input: CpG report filtered by specific chromosome
# output: scatter plots of + and - strand that chromosome
# steps: 
# 1) read the file 
# 2) keep track of strand (+ or -) and position
# 3) scatter plot of + positions in one color and - positions in another
#       across the chromosome

import os
import matplotlib.pyplot as plt

# input parameters
#inputfile = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0055_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.chrX.txt"
#inputfile = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0083_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.chrX.txt"
inputfile = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0088_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.chrX.txt"
#inputfile = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0095_Placenta_Sample1_trimmed_R1_bismark_bt2_pe.deduplicated.chrX.txt"
outputpng = inputfile.replace(".txt",".png")
outputtxt = inputfile.replace(".txt","strand.txt")

# load data from file from chrX file
print ("Loading input data")
inputhandle = open(inputfile,"r")
inputdata = inputhandle.readlines();
inputhandle.close()

# hard coded positions of data we want to plot
strand_index = 1
position_index = 3

# list of positions that have + or - strand
plus_strand_events = []
minus_strand_events = []

for input in inputdata:
    input = input.replace("\n","")
    items = input.split("\t")
    if (items[strand_index] == "+"):
        plus_strand_events.append(int(items[position_index]))
    elif (items[strand_index] == "-"):
        minus_strand_events.append(int(items[position_index]))
    if ((len(plus_strand_events) + len(minus_strand_events)) % 1000000 == 0):
        print ("Processing...")

# output + and - strand count in text file
print ("Writing output file")
outputhandle = open(outputtxt,"w")
print("\t".join(["Number Plus Strand","Number Minus Strand", "Fold Change Num_minus/Num_plus"]),file = outputhandle)
num_plus = len(plus_strand_events)
num_minus = len(minus_strand_events)
fc = num_minus / num_plus
print("\t".join([str(num_plus),str(num_minus), str(fc)]),file = outputhandle)
outputhandle.close()
print ("Created output: ", outputtxt)

# plot + and - strand count
print ("Making figure")

fig, ax = plt.subplots(figsize=(4, 3))

# x coordinates of scatter plot will be positions
# y coordinates of plus strand events will be 2
# y coordinates of minus strand events will be 1
# plus_strand_ycoordinates = [2] * len(plus_strand_events)
# minus_strand_ycoordinates = [1] * len(minus_strand_events)
#ax.scatter(x = plus_strand_events, y = plus_strand_ycoordinates, marker = '.', s = 6, c = "red")
#ax.scatter(x = minus_strand_events, y = minus_strand_ycoordinates, marker = '.', s = 6, c = "blue")

bar_labels = ["Number Plus Strand", "Number Minus Strand"]
bar_height = [len(plus_strand_events),len(minus_strand_events)]
ax.bar(x = bar_labels, height = bar_height, width = 0.5)

plt.savefig(outputpng, bbox_inches='tight')
print ("Created plot: ", outputpng)
