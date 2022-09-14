# run with module python/3.7.1 so the order is maintained

import json

data = json.load(open("TrimmedFastqFilePaths_local.rdgroups.json","r"))

germline_samples = data["all_germline_samples"]

for sample in data.keys():
    if ("trimmed" in sample and not sample in germline_samples):
        sample_parts = sample.split("_")
        sampleid = sample_parts[0]
        if (sampleid[-1].isalpha()):
            sampleid = sampleid[0:len(sampleid)-1]
        find_normal = [i for i in germline_samples if sampleid+"Blood" in i or sampleid+"BL" in i]
        if (len(find_normal) > 0):
            data[sample]["normal"] = find_normal[0]

new_data = json.dumps(data,indent = 4)

outputhandle = open("TrimmedFastqFilePaths_local.rdgroups.withnormals.json","w")

print (new_data, file = outputhandle)

