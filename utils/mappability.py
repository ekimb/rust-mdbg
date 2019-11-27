import sys
import parse_sequences_file

k, l, kmer_to_seq, kmer_abundance = parse_sequences_file.parse(sys.argv[1])

#print(k,l)

nb_unique = len([kmer for kmer in kmer_abundance if kmer_abundance[kmer] == 1])
nb_distinct = len(kmer_abundance)
perc_only_once = 100.0*nb_unique / nb_distinct

print("percentage of distinct %d-mers seen only once: %.2f" % (k,perc_only_once))

if perc_only_once == 100.0:
    exit(1)
