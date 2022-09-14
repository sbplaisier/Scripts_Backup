# compare fastas of gametologs

# import packages
from Bio import pairwise2
from Bio import Align

# load gametologs
gametologs = { "DDX3X":"DDX3Y",
    "PCDH11X":"PCDH11Y",
    "USP9X":"USP9Y",
    "ZFX":"ZFY",
    "UTX":"UTY",
    "XIST":"SRY",
    "KDM5C":"KDM5D",
    "PRKX":"PRKY",
    "RPS4X":"RPS4Y1",
    "EIF1AX":"EIF1AY",
    "NLGN4X":"NLGN4Y",
    "TGIF2LX":"TGIF2LY",
    "SOX3":"SRY",
    "AMELX":"AMELY",
    "TBL1X":"TBL1Y",
    "DBX":"DBY",
    "TMSB4X":"TMSB4Y",
    "CXorf15":"CYorf15A",
    "CXorf15":"CYorf15B",
    "SMCX":"SMCY",
    "RPS4X":"RPS4Y2",
    "VCX":"VCY",
    "RBMX":"RBMY" }

# load fastas
#chrXfasta = "gametology_chrX_genes_compare.fa"
#chrYfasta = "gametology_chrY_genes_compare.fa"
chrXfasta = "gametology_chrX_exons_compare_fused.fa"
chrYfasta = "gametology_chrY_exons_compare_fused.fa"

chrXhandle = open(chrXfasta,"r")
chrXseq = chrXhandle.readlines()
chrXhandle.close()

chrYhandle = open(chrYfasta,"r")
chrYseq = chrYhandle.readlines()
chrYhandle.close()

# set output file
outputfile = "gametology_length_similarity.csv"
outputhandle = open(outputfile,"w")
#print("Pair,length,length_difference,similarity",file=outputhandle)
print("Pair,length,length_difference",file=outputhandle)

# iterate through them
for chrXg in gametologs:
    chrXseqfound = ""
    chrYseqfound = ""
    chrYg = gametologs[chrXg]
    for chrX_index in range(0,len(chrXseq)-1):
        chrXline = chrXseq[chrX_index]
        if (">" in chrXline and chrXg in chrXline):
            chrXseqfound = chrXseq[chrX_index + 1]
        for chrY_index in range(0,len(chrYseq)-1):
            chrYline = chrYseq[chrY_index]
            if (">" in chrYline and gametologs[chrXg] in chrYline):
                chrYseqfound = chrYseq[chrY_index + 1]
    if (chrXseqfound != "" and chrYseqfound != ""):
        gametolog_pair = chrXg + ":" + gametologs[chrXg]
             
        # print gene length for each pair
        length_pair = str(len(chrXseqfound)) + ":" + str(len(chrYseqfound))
        length_diff = len(chrXseqfound) - len(chrYseqfound)

        print (gametolog_pair + "," + length_pair + "," + str(length_diff), file = outputhandle)

        # print out sequences to submit for alignment online
        outseqX = open("align_exons_" + chrXg + ".txt","w")
        print (chrXseqfound,file=outseqX)
        outseqX.close()

        outseqY = open("align_exons_" + gametologs[chrXg] + ".txt","w")
        print (chrYseqfound,file=outseqY)
        outseqY.close()
        
'''
        # print sequence similarity for each pair
        percent_match = 0

        match = 2
        mismatch = -1
        gap_open = -10
        gap_extend = -0.5

        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend
        aligner.mode= 'local'
        for alignment in aligner.align(chrXseqfound, chrYseqfound):
            print("Score = %.1f:" % alignment.score)
        print (gametolog_pair + "," + length_pair + "," + str(length_diff) + "," + str(percent_match), file = outputhandle)
'''

outputhandle.close()
