""" purpose: 
    report which % of each gene is found in the graph
    hits format: [hit_id] [hit_pos] [list of (gene_id,gene_pos)]
    reference file: many lines, each line is: [seq_id] [list of minimizers]
"""

import sys
gene_filename = sys.argv[1]
hits_filename = sys.argv[2]

k=10 

from collections import defaultdict

gene_pos_covered = dict()
for line in open(gene_filename):
    ls = line.split()
    gene_id = ls[0]
    minimizers = tuple(map(int,ls[1:]))
    if len(minimizers) < k: continue
    nb_minimizers = len(minimizers)-k+1
    gene_pos_covered[gene_id] = [False]*nb_minimizers
    if nb_minimizers <= 0:
        print("uh oh",nb_minimizers,gene_pos_covered[gene_id])


for line in open(hits_filename):
    ls = line.split()
    hit_id, hit_pos = ls[:2]
    gene_list = eval(" ".join(ls[2:]))
    for gene_id, gene_pos in gene_list:
        gene_pos_covered[gene_id][gene_pos]=True


for gene_id in gene_pos_covered:
    nb_covered=sum(gene_pos_covered[gene_id])
    nb_mins=len(gene_pos_covered[gene_id])
    print(gene_id,nb_mins,"%.2f" % (nb_covered/nb_mins*100.0))
