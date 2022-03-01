# read output of samtools depth
# for each readdepth.txt file
# keep sum and running tally
# divide sum/count for average at the end
# write out results with .chr extension

import os

inputdirectory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/avg_read_depth/rna/"
files = os.listdir(inputdirectory)
chrrdfiles = [f for f in files if ("readdepth.chr.txt" in f and "Placenta" in f)]
rdfiles = [f for f in files if ("readdepth.txt" in f and "Placenta" in f and not f.replace(".txt",".chr.txt") in chrrdfiles)]

for rdfile in rdfiles:
#rdfile = rdfiles[0]
    print ("Processing " + inputdirectory + rdfile)
    fullrdfile = inputdirectory + rdfile

    depthcount = {}
    counter = 0
    with open(fullrdfile) as f:
        line = f.readline()

        while line: 
            counter += 1

            line = line.replace("\n","")
            (chrom,pos,depth) = line.split("\t")

            if (not chrom in depthcount):
                depthcount[chrom] = (int(depth),1)
                print ("Starting " + chrom)
            else:
                (current_sum, current_count) = depthcount[chrom]
                current_sum += int(depth)
                current_count += 1
                depthcount[chrom] = (current_sum,current_count)

            #if (counter == 5):
            #    break
            if (counter % 10000000 == 0):
                print ("At position: " + str(counter))

            line = f.readline()

    outputfile = fullrdfile.replace(".txt",".chr.txt")
    outputhandle = open(outputfile, "w")
    for c in depthcount:
        (current_sum,current_count) = depthcount[c]
        print (c + "\t" + str(current_sum/current_count), file = outputhandle)
    outputhandle.close()
