# read output of samtools depth
# figure out of max of column 3 (read depth for each gene)
# figure out the top 10

import os

inputdirectory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/avg_read_depth/"
inputfile = "gencode.v29.annotation.gene.maxreaddepth.OBG0055_Plac1b.txt"

print ("Processing " + inputdirectory + inputfile)

top10 = []
temp = "start\tstart\t0"
top10.append(temp)

counter = 0
with open(inputdirectory+inputfile) as f:
    line = f.readline()

    while line: 
        counter += 1

        line = line.replace("\n","")
        linesplit = line.split("\t")
        
        currentmaxline = top10[0]
        currentmaxlinesplit = currentmaxline.split("\t")
        currentmax = currentmaxlinesplit[2]
        linemax = linesplit[2]
        if ("b" in linemax):
            linemax = 0
        else:
            linemax = int(linemax)
        
        if (linemax > int(currentmax)):
            top10.insert(0, line)
        if (len(top10) > 10):
            top10 = top10[0:10]
            
        #if (counter == 30):
        #    break
        if (counter % 10000000 == 0):
            print ("At position: " + str(counter))

        line = f.readline()

outputfile = inputfile.replace(".txt",".max.txt")
outputhandle = open(inputdirectory+outputfile, "w")
print ("\n".join(top10), file=outputhandle)
outputhandle.close()
