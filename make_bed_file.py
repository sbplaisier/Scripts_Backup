# write out bed file format for every entry with filter term in the filter field 

import os
import sys

annotation_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/avg_read_depth/gencode.v29.annotation.exon.txt"
#annotation_file = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/avg_read_depth/gencode.v29.annotation.gene.txt"
#annotation_file = sys.argv[1] 
indices_to_print = [0,3,4,8]
#indices_to_print = [2,3,3,0]
filter_term = "exon"
#filter_term = "gene"
#filter_term = "Z"
filter_index = 2
#filter_index = 4

outputfile = annotation_file.replace(".txt",".bed")
outputhandle = open(outputfile, "w")

counter = 0
with open(annotation_file) as f:
    line = f.readline()

    while line: 
        counter += 1

        if ("chr" in line):

            line = line.replace("\n","")
            linesplit = line.split("\t")
            #print (linesplit)

            print_string = ""
            if (filter_term in linesplit[filter_index]):
                subset = [linesplit[i] for i in indices_to_print]
                print_string = "\t".join(subset)
                print(print_string, file=outputhandle)

            #if (counter == 100000):
            #    break
            #if (counter % 10000000 == 0):
            #    print ("At position: " + str(counter))

        line = f.readline()

outputhandle.close()
