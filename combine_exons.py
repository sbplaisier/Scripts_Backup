# add gene name and exon number to header of fasta

# load fastas and matching gtf file
fasta_file = "gametology_chrX_exons_compare_annotated.fa"
#fasta_file = "gametology_chrY_exons_compare_annotated.fa"

inputhandle = open(fasta_file,"r")
fasta_data = inputhandle.readlines()
inputhandle.close()

# combine the sequence for multiple exons in the same gene
#  make sure to check for repeats

# hard-coded
range_index = 0
gene_name_index = 2
exon_index = 3

# use hash to store gene sequence
# key = gene_name
# value = list of exon numbers, list of sequence
collect_exons = {}
for i in range(0,len(fasta_data)-1):
    line = fasta_data[i]
    line = line.replace("\n","")
    if (">" in line): # header
        # parse header
        (coordinates, chromosome, gene_name, exon_number) = line.split(" ")
        gene_name = gene_name.replace("gene_name:","")
        exon_number = exon_number.replace("exon_number:","")        
        sequence = fasta_data[i+1]
        sequence = sequence.replace("\n","")
        i += 1

        if(not gene_name in collect_exons):
            exons = {}
            exons["exon_number"] = [int(exon_number)]
            exons["seq"] = [sequence]
            collect_exons[gene_name] = exons
        else: 
            exons = collect_exons[gene_name]
            if (not int(exon_number) in exons["exon_number"]):
                exons["exon_number"].append(int(exon_number))
                exons["seq"].append(sequence)
            collect_exons[gene_name] = exons

# iterate through all entries stringing together the sequences in order
combined_exons = {}
for gene in collect_exons:
    original_exon_array = collect_exons[gene]["exon_number"]
    sorted_exon_array = original_exon_array.copy()
    sorted_exon_array.sort()

    full_sequence = ""
    exons_added = []

    # find the index in the original exon_number array
    # tack on the sequence at that index

    for i in sorted_exon_array:
        exon_original_index = original_exon_array.index(i)
        exons_added.append(i)
        full_sequence += collect_exons[gene]["seq"][exon_original_index]

    combined_exons[gene] = (exons_added,full_sequence)

# set output file
outputfile = fasta_file.replace("_annotated.fa","_fused.fa") 
outputhandle = open(outputfile,"w")

# print all combined sequences as fasta

for gene in combined_exons:
    (exon_list,fullseq) = combined_exons[gene]

    print ("> " + gene + " exons:" + str(exon_list), file = outputhandle)
    print (fullseq, file = outputhandle)


outputhandle.close()
