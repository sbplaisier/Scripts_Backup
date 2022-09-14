import sys
import re

inputfile = sys.argv[1]

handle = open(inputfile,"r")
data = handle.readlines()
handle.close()


map_qualities = []

for line in data:
    m = re.search("MQ:\w:\d",line)
    if (m):
        mapq_string = m.group()
        mapq_split = mapq_string.split(":")
        mq = int(mapq_split[-1])
        map_qualities.append(mq)

outputfile = inputfile.replace(".txt","_MQ.txt")
outhandle = open(outputfile,"w")

max_mq = max(map_qualities)
min_mq = min(map_qualities)
print("Max map quality: " + str(max_mq),file = outhandle)
print("Min map quality: " + str(min_mq),file = outhandle)

print("Breakdown: ",file = outhandle)
for i in range(min_mq,max_mq):
    counter = 0
    for m in map_qualities:
        if (m == i):
            counter += 1
    print("\t"+str(counter) + " entries had " + str(i) + " quality", file = outhandle)

outhandle.close()
