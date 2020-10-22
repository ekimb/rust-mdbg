import sys
if len(sys.argv) < 3 or ".sequences" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: [graph.sequences] [final.sequences] \n will take info from [graph.sequences] (that should contain only kmers) and fill in [final.sequences]\n note: overwrites the third column of [final.sequences]")

# read [origin.sequences] file
kmer_to_seq = dict()
k, l = 0, 0
for line in open(sys.argv[1]):
    # ignore #'s lines except for getting the k value
    if line.startswith('#'):
        if line.startswith('# k = '):
            k = int(line.split()[-1])
        if line.startswith('# l = '):
            l = int(line.split()[-1])
        continue
    spl = line.split()
    id = spl[0]
    minims = tuple(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-5]))
    seq = spl[-5]
    kmer_to_seq[minims] = seq

# read and cache [final.sequences]
final_sequences_file = []
for line in open(sys.argv[2]):
    final_sequences_file += [line]

import string
old_chars = "ACGT"
replace_chars = "TGCA"
tab = str.maketrans(old_chars,replace_chars)
def revcomp(s):
    return s.translate(tab)[::-1]

from itertools import zip_longest # for Python 3.x
def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def double_every_k(k, iterable):
    # weird double so that we get kmer that share the same minimizer as the last
    # double_every_k(3, abcdefg) -> (abc,cde,efg)
    # to test: print(list(double_every_k(5,list(range(20)))))
    counter = 1
    for elt in iterable:
        if counter > 0 and counter % k == 0:
            counter = 1
            yield elt
        yield elt
        counter += 1

for i in range(len(final_sequences_file)):
    line = final_sequences_file[i]
    if line.startswith('#'): continue
    spl = line.split()
    utg = spl[0]
    spl_end = -3
    minims = tuple(map(lambda x: int(x.strip('(').strip(')').replace(',','')),spl[1:spl_end]))
    seq = spl[spl_end]
    
    #magic happens here
    
    # don't wanna do it k-mer-per-k-mer, let's be more rough than that
    #for j in range(len(minims)-k+1):
        #kmer = minims[j:j+k]
    whole_seq = ""
    for kmer in grouper(k, double_every_k(k, minims)):
        do_revcomp = False
        if None in kmer: continue 
        if kmer not in kmer_to_seq:
            kmer = kmer[::-1]
            do_revcomp = True
        if kmer not in kmer_to_seq:
            print("kmer not found",kmer)
            exit(1)
        seq = kmer_to_seq[kmer]
        if do_revcomp:
            seq = revcomp(seq)
        #print(seq)
        if len(whole_seq) == 0:
            whole_seq = seq
        else:
            if whole_seq[-l:] != seq[:l]:
                """
                # what might happen is when a repetition of the same minimizer occurs
                # whole_seq is, in terms of minimizers: a,b,c,d,X1,X2
                # seq is X1,e,f,g,h
                # but for some reason, seq has a sequencing error and X2 is missed
                # then gotta patch dirtily
                print("problem with overlap") #todo understand what's up
                print(whole_seq[-l:])
                print(seq[:l])
                print(whole_seq)
                print(seq)
                print(minims)
                exit(1)
                """
            whole_seq += seq[l:]

    final_sequences_file[i] = "%s\t%s\t%s\tPLACEHOLDER\tPLACEHOLDER\n" % (utg, list(minims), whole_seq)

# write [final.sequences] with sequence info 
output = open(sys.argv[2],'w')
for line in final_sequences_file:
    output.write(line)
output.close()
