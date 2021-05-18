""" purpose: scans stdin for kminmers found in a reference file
    stdin format: [seq_id] [list of minimizers]
    reference file: many lines, each line is: [seq_id] [list of minimizers to extract kminmers from]
"""

import sys
genome_minspace_filename = sys.argv[1]
k=10 

graph_mode = "-g" in sys.argv

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
        kminmers[kminmer]       += [(seq_id,i)]
        kminmers[kminmer[::-1]] += [(seq_id,i)]
        #kminmers_sets.add(set(kminmer))
        assert(len(kminmer)==10)
        # hack for speed
        kminmer_str     = str(list(kminmer))
        kminmer_str_inv = str(list(kminmer[::-1]))
        #print(kminmer_str,kminmer_str_inv)
        kminmers[kminmer_str]     += [(seq_id,i)]
        kminmers[kminmer_str_inv] += [(seq_id,i)]

if graph_mode:
    for line in sys.stdin:
        if line[0] == "#": continue
        # 16606   [27472887960080780, 26945328166221359, 83024137861838436, 183455804785478733, 54911163836344167, 170342695321694208, 91112118779713090, 54911163836344167, 83024137861838436, 183455804785478733]       GCCGAGAGGCTGAAGGCGCTCCCCTGCTAAGGGAGTATGCGGTCAA        AAGCTGCATCCGGGGTTCGAATCCCCGCCTCACCGCCATTTGCATCCGTAGCTCAGCTGGATAGAGTACTCGGCTACGAACCGAGCGGTCGGAGGTTCGAATCCTCCCGGATGCACCATATTCTACGTACTTTCAGCGATGAAGGTATGGAAGAGGTGGCGGTATAACCGCAGGCACCAGGGAGGATAACGTTGCTTTAGCAACGGCCCGAAGGGCGAGCCGCAAGGCGAGTAATCCTCCCGGATGCACCATCT        CTTACTTGATATGGCTTTAGTAGCGATATCAATATCAGCAGTAAAATAAATTTTCCCGATGCATCCGTAGCTCAGCTGGATAGAGTACTCGGCTACGAACCGAGCGGTC   *       *       (24, 33)
        ls = line.split()
        #kminmer = tuple(eval(" ".join(ls[1:k+1])))
        kminmer = " ".join(ls[1:k+1])
        if kminmer in kminmers:
            print("*","*",kminmers[kminmer] + kminmers[kminmer[::-1]])
else:
    for line in sys.stdin:
        if ":" in line or "Parsing" in line: continue # debug stuff output by mdbg
        ls = line.split()
        seq_id, minimizers = ls[0], tuple(map(int,ls[1:]))
        for i in range(len(minimizers)-k+1):
            kminmer = minimizers[i:i+k]
            if kminmer in kminmers or kminmer[::-1] in kminmers:
                #print("hit!",seq_id,"at pos",i,"found in ref(s)",kminmers[kminmer] + kminmers[kminmer[::-1]])
                print(seq_id,i,kminmers[kminmer] + kminmers[kminmer[::-1]])
