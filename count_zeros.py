# input: 
#   2) directory of allele balance info
#  output: tsv with another column added for call rate
#  steps:
#   1) if column header contains 'count'
#   2) look through that column and report how many items in that column are 0

import os

#ase_directory = "/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/all_samples_joint/asereadcounter/phased_allele_balance/"
ase_directory = "/data/CEM/wilsonlab/projects/placenta/PlacentaXCI_XXfemales_firstbatch/04_phasing/phased_allele_balance/"

allele_balance_files = list()
for (dirpath, dirnames, filenames) in os.walk(ase_directory):
    allele_balance_files += [os.path.join(dirpath, file) for file in filenames if "allele_balance.tsv" in file or "allele_balance_data.tsv" in file]

if (len(allele_balance_files) == 0):
    print ("Didn't find any allele balance files here")

# check to see what the delimiter of the files is
# sometimes file is labeled tsv when it is really a csv
inputhandle = open (allele_balance_files[0], "r")
lines = inputhandle.readlines()
inputhandle.close()
headerline = lines[0]
headers_tab = headerline.split("\t")
headers_comma = headerline.split(",")
delimiter = "\t"
if (len(headers_comma) > len(headers_tab)):
    delimiter = ","

for allele_balance_file in allele_balance_files:
    inputhandle = open (allele_balance_file, "r")
    lines = inputhandle.readlines()
    inputhandle.close()
    
    outputfile = allele_balance_file.replace(".tsv","_count0.tsv")
    outputhandle = open (outputfile, "w")

    headerline = lines[0]
    headers = headerline.split(delimiter)
    count_indices = []
    count_headers = []
    for i in range(0,len(headers)):
        if ('count' in headers[i]):
            count_indices.append(i)
            count_headers.append(headers[i])
        
    print ("\t".join(count_headers), file=outputhandle)
    
    counts = [0] * len(count_headers) # keep track of the counts

    for line in lines[1:len(lines)]:
        items = line.split(delimiter)
        for i in range(0,len(count_indices)):
            if (items[count_indices[i]] == "0"):
                counts[i] += 1

    counts_str = [str(x) for x in counts]
    percents_str = ['{:.2f}%'.format( 100 * x / (len(lines) - 1)) for x in counts]
    print ("\t".join(counts_str),file=outputhandle)
    print ("\t".join(percents_str),file=outputhandle)
    outputhandle.close()
