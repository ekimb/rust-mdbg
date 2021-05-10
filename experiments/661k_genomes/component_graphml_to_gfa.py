# it's a bit cheating here : uses the original gfa the graphml component is from

import gzip
orig_gfa = "../reggraph-k10-p0.001-l12.u.gfa.gz"


import sys
graphml_file = sys.argv[1]

unitigs = set()
for line in open(graphml_file):
    if "node id" in line:
        ls = line.split()
        node = ls[1].strip('id=\"').strip('\">')
        #print(node)
        unitigs.add(node)

for line in gzip.open(orig_gfa):
    line = line.decode().strip()
    ls = line.split()
    if line.startswith('S'):
        if ls[1] in unitigs:
            print(line)
    elif line.startswith('A'):
        if ls[1] in unitigs:
            print(line)
    elif line.startswith('L'):
        if ls[1] in unitigs and ls[3] in unitigs:
            print(line)
    else:
        print(line)
