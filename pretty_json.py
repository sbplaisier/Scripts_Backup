import json
import sys

inputfile = sys.argv[1]

with open(inputfile, "r") as read_file:
    data = json.load(read_file)

outputfile = inputfile.replace(".json", "_formatted.json")
print(json.dumps(data, indent=4, sort_keys=True), file=open(outputfile,"w"))
