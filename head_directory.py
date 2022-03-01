import os

directory = "/data/CEM/wilsonlab/lab_generated/placenta/2021_placenta-decidua_testing/WGBS"
outputfile = "/home/splaisie/Placenta_Data_Processing/WGBS_heads.txt"

listOfFiles = list()
for (dirpath, dirnames, filenames) in os.walk(directory):
    listOfFiles += [os.path.join(dirpath, file) for file in filenames if ".gz" in file]

for file in listOfFiles:
    echocommand = "echo " + file + " >> " + outputfile
    print (echocommand)
    systemcommand = "gzip -cd " + file + " | head -n 1 >> " + outputfile
    print (systemcommand)
    #os.system(systemcommand)
    
