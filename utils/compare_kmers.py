import sys
if len(sys.argv) < 3 or ".sequences" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: <file1.sequences> <file2.sequences> [genome.fasta]\n will compare the set of kmers in those two files\nand optionnally analyze the missing kmers across the refgenome\n")

# tip: make file1 be the .sequences file constructed from a reference genome
# and file2 be one made from the reads
file1= sys.argv[1]
file2= sys.argv[2]

from parse_sequences_file import parse
k, l, node_minims,  kmer_seq1, osef, origins1 = parse(file1)
osef, osef, node_minims2, kmer_seq2, osef, origins2 = parse(file2)

kmers1 = set(kmer_seq1.keys())
kmers2 = set(kmer_seq2.keys())

print(len(kmers1),"kmers in",file1)
print(len(kmers2),"kmers in",file2)

#print("intersection:",len(kmers1 & kmers2))

kmers1_in_kmers2 = len([x for x in kmers1 if x in kmers2])
kmers1_not_in_kmers2  = len(kmers1) - kmers1_in_kmers2
print("number of kmers from %s that are in %s" % (file1,file2), kmers1_in_kmers2, "(%.2f)%%" % ((1.0*kmers1_in_kmers2)/len(kmers1)*100.0),",", kmers1_not_in_kmers2,"are not")

# use genome to find kmer location
debug_missing_kmers = len(sys.argv) > 3
if debug_missing_kmers:
    
    """
    # unused code
    def rev_compl(st):
        nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(nn[n] for n in reversed(st))
    from pyfasta import Fasta
    genome_filename = sys.argv[3]
    g = Fasta(genome_filename)
    genome = g[list(g.keys())[0]] # assume a single chromosome
    genome_rc = rev_compl(genome)
    assert(len(g.keys()) == 1)
    """

    print("the",kmers1_not_in_kmers2,"missing kmers:")
    for kmer in [x for x in kmers1 if x not in kmers2]:
        if kmer in origins1:
            origin = origins1[kmer]
            f = open("missing_kmers.txt", "a")
            print(kmer,origin)
        else:
            print(kmer,"not found in reference '.sequences' file, but should have been")
        if kmer[::-1] in origins1:
            print("Reverse complement in reference")
    f.close()
