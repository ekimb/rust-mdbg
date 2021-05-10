""" purpose: 
    find which graph node corresponds to which genome
"""

import sys
hits_filename = "hits.txt"
utg_kinmer_filename = "all_kminmers.unitig_assoc.txt"

k=10 

from collections import defaultdict

kminmer_dict = defaultdict(list)
unitigs_dict = defaultdict(set)
for line in open(hits_filename):
    ls = line.split()
    hit_id, hit_pos = ls[:2]
    kminmer_list = eval(" ".join(ls[2:]))
    for kminmer_and_pos in kminmer_list:
        kminmer = kminmer_and_pos[0]
        kminmer_dict[kminmer] += [hit_id]

for line in open(utg_kinmer_filename):
    unitig, kminmer = line.split()
    unitigs_dict[unitig] |= set(kminmer_dict[kminmer])


nb_genomes_per_unitig = 0 
for unitig in unitigs_dict:
    nb_genomes_per_unitig += len(unitigs_dict[unitig])

nb_genomes_per_unitig /= len(unitigs_dict)
#print(nb_genomes_per_unitig)

for unitig in unitigs_dict:
    print(unitig,unitigs_dict[unitig])
