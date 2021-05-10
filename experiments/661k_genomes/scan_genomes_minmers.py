""" purpose: scans stdin for kminmers found in a reference file
    stdin format: [seq_id] [list of minimizers]
    reference file: many lines, each line is: [seq_id] [list of minimizers to extract kminmers from]
"""

import sys
genome_minspace_filename = sys.argv[1]
k=10 

from collections import defaultdict
kminmers = defaultdict(list)
#kminmers_sets = set()

for line in open(genome_minspace_filename):
    line = line.replace('[','').replace(']','').replace(',','')
    ls = line.split()
    seq_id = ls[0]
    minimizers = tuple(map(int,ls[1:]))
    if len(minimizers) < k: continue
    for i in range(len(minimizers)-k+1):
        kminmer = minimizers[i:i+k]
        kminmers[kminmer] += [(seq_id,i)]
        #kminmers_sets.add(set(kminmer))
        assert(len(kminmer)==10)

for line in sys.stdin:
    if ":" in line or "Parsing" in line: continue # debug stuff output by mdbg
    ls = line.split()
    seq_id, minimizers = ls[0], tuple(map(int,ls[1:]))
    for i in range(len(minimizers)-k+1):
        kminmer = minimizers[i:i+k]
        if kminmer in kminmers or kminmer[::-1] in kminmers:
            #print("hit!",seq_id,"at pos",i,"found in ref(s)",kminmers[kminmer] + kminmers[kminmer[::-1]])
            print(seq_id,i,kminmers[kminmer] + kminmers[kminmer[::-1]])
