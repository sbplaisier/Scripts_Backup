# add gene name and exon number to header of fasta

# load fastas and matching gtf file
#fasta_file = "gametology_chrX_exons_compare.fa"
#gtf_file = "gametology_chrX_exons_compare.gtf"
fasta_file = "gametology_chrY_exons_compare.fa"
gtf_file = "gametology_chrY_exons_compare.gtf"

inputhandle = open(fasta_file,"r")
fasta_data = inputhandle.readlines()
inputhandle.close()

inputhandle = open(gtf_file,"r")
gtf_data = inputhandle.readlines()
inputhandle.close()

# iterate through lines of gtf file filling a hash, annotation
#  key = start-end, value = [chromosome,gene_name,exon_number]

# hard-coded
chr_index = 0
exon_index = 2
start_index = 3
end_index = 4
info_index = 8

annotation = {}
for gtf_line in gtf_data:
    gtf_elements = gtf_line.split("\t")
    chromosome = gtf_elements[chr_index] 
    exon_check = gtf_elements[exon_index]
    start = gtf_elements[start_index]
    end = gtf_elements[end_index]
    info = gtf_elements[info_index]
    info_elements = info.split(";")
    gene_name = ""
    exon_number = ""

    for ie in info_elements: 
        if ("gene_name" in ie):
            gene_name = ie.replace("gene_name \"","")
            gene_name = gene_name.replace("\"","")
            gene_name = gene_name.replace(" ","")
        if ("exon_number" in ie):
            exon_number = ie.replace("exon_number ","")
            exon_number = exon_number.replace("\"","")
            exon_number = exon_number.replace(" ","")
            
    annotation[start+'_'+end] = [chromosome,"gene_name:"+gene_name,"exon_number:"+exon_number]

# set output file
outputfile = fasta_file.replace(".fa","_annotated.fa") 
outputhandle = open(outputfile,"w")

# read fasta file and print with annotation info in header
for i in range(0,len(fasta_data)-1):
    line = fasta_data[i]
    line = line.replace("\n","")
    if (">" in line): # header
        # parse header chr_:start-end
        (chr_found,range_found) = line.split(":")
        (start,end) = range_found.split("-")
        sequence = fasta_data[i+1]
        sequence = sequence.replace("\n","")
        i += 1

        line_annotated = line + " " + " ".join(annotation[str(int(start)+1)+"_"+end]) # has different base so start in fasta is one less than start in gtf

        print(line_annotated, file=outputhandle)
        print(sequence, file=outputhandle)

outputhandle.close()
