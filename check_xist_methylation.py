# check all CpG methylation extraction for sites in XIST (chrX:73,820,656-73,852,723)

meth_files = ["/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0023_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0023_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0055_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0083_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0085_Placenta_Sample1_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0088_Placenta_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/methylation3/CpG_context_OBG0095_Placenta_Sample1_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/CpG_context_OBG0023_Decidua_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/CpG_context_OBG0055_Decidua_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/CpG_context_OBG0083_Decidua_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/CpG_context_OBG0085_Decidua_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt",
"/data/CEM/wilsonlab/projects/placenta/2021_placenta_test/decidua/wgbs/CpG_context_OBG0088_Decidua_Merge_trimmed_R1_bismark_bt2_pe.deduplicated.txt"]

#A00124:297:H7V5VDSX2:2:1101:6424:1000_1:N:0:TTCAATAG+TCGTGGGA   +       chr2    46220329        Z

xist_chr = "chrX"
xist_start = 73820656
xist_end = 73852723

for meth_file in meth_files:
    print (meth_file)

    output_file = meth_file.replace(".txt","_xist.txt")
    outhandle = open(output_file,"w")
    
    counter = 0
    unmeth_counter = 0
    meth_counter = 0
    with open(meth_file) as f:
        line = f.readline()

        while line: 
            counter += 1

            line = line.replace("\n","")
 
            if ("chr" in line):
                line_split = line.split("\t")
                chr_found = line_split[2]
                pos_found = line_split[3]
                meth_found = line_split[4]
                if (chr_found == xist_chr and int(pos_found) >= xist_start and int(pos_found) <= xist_end):
                    print(line,file=outhandle)
                    if (meth_found == 'z'):
                        unmeth_counter += 1
                    elif(meth_found == 'Z'):
                        meth_counter += 1


            line = f.readline()
            
    print ("\nSummary: ",file=outhandle)
    print ("unmethylated Cs: " + str(unmeth_counter),file=outhandle)
    print ("methylated Cs: " + str(meth_counter),file=outhandle)
    print ("percent methylated xist = " + str((meth_counter)/(meth_counter + unmeth_counter)),file=outhandle)
    outhandle.close() 

