import sys

ase_file = sys.argv[1]
vcf_file  = sys.argv[2]

inputhandle = open(ase_file,"r")
lines = inputhandle.readlines()
inputhandle.close()
ase_positions =[]
for line in lines[1:len(lines)]:
    items = line.split("\t")
    ase_positions.append(items[1])

vcf_position_index = 1
if (".vcf" in vcf_file):
    vcf_position_index = 1
if (".bed" in vcf_file):
    vcf_position_index = 2

inputhandle = open(vcf_file,"r")
lines = inputhandle.readlines()
inputhandle.close()
store_filtered_variants =[]
for line in lines:
    if (not "#" in line):
        items = line.split("\t")
        store_filtered_variants.append(items[vcf_position_index])

matches = 0
for ase_pos in ase_positions:
    if ase_pos in store_filtered_variants:
        matches += 1

print (str(matches/len(ase_positions) * 100) + "%")
