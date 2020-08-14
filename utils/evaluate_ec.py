"""
 this script takes as input:
 - the .ec_data file for a reference genome (processed by rust-mhdbg with min-abundance of 1)
 - the .ec_data file for a set of reads 
 - and optionnally,
   the .ec_data file for the same set of reads, processed differently (e.g. corrected)
 - and also optionnally,
   the .poa.ec_data file that contains which reads are retrieved for the template
 and outputs the needleman-wunch alignment of each read to the reference (semi-global aln)  
 and optionally outputs a comparison between the two set of reads (e.g. corrected vs uncorrected)
"""

nb_processes = 8

import sys
if len(sys.argv) < 3 or ".ec_data" not in sys.argv[2] or ".ec_data" not in sys.argv[1]:
    exit("input: <reference.ec_data> [reads.ec_data] [reads.corrected.ec_data] [reads.poa.ec_data]\n will evaluate accuracy of minimizers in reads\n")

import evaluate_poa

file1 = sys.argv[1]
file2 = sys.argv[2]
if len(sys.argv) > 3:
    file3 = sys.argv[3]
if len(sys.argv) > 4:
    file_poa = sys.argv[4]
    poa_d, poa_d_itv, poa_reads = evaluate_poa.prepare_eval_poa(file_poa)
else:
    file_poa = None

def parse_file(filename):
    res = []
    counter = 0
    seq_id = ""
    for line in open(filename):
        if counter % 5 == 0:
            seq_id = line.strip()
        if counter % 5 != 2:
            counter += 1 
            continue
        spl = line.split()
        minims = list(map(int,spl))
        res += [(seq_id,minims)]
        counter += 1
    return res


reference = parse_file(file1)
assert(len(reference)==1)
reads  = parse_file(file2)
reads2 = parse_file(file3) if len(sys.argv) > 3 else None

print("loaded",len(reference),"reference,",len(reads),"reads")

reference = reference[0][1]

# adapted from NW code here https://stackoverflow.com/questions/2718809/how-to-diff-align-python-lists-using-arbitrary-matching-function
# 
def semiglobal_align(a, b):
    # This is regular Needleman-Wunsch scoring.
    # which is not quite edit distance: match would need to be 0 for ED.
    # But i dunno how to compute semiglobal edit distance so i'll take that for now!
    replace_func = lambda x,y: 1 if x==y else -1
    insert = -1
    delete = -1

    #for traceback
    ZERO, LEFT, UP, DIAGONAL = 0, 1, 2, 3

    len_a = len(a)
    len_b = len(b)

    matrix = [[(0, ZERO) for x in range(len_b + 1)] for y in range(len_a + 1)]

    """
    #this is for NW
    for i in range(len_a + 1):
        matrix[i][0] = (insert * i, UP)

    for j in range(len_b + 1):
        matrix[0][j] = (delete * j, LEFT)
    """

    for i in range(1, len_a + 1):
        for j in range(1, len_b + 1):
            replace = replace_func(a[i - 1], b[j - 1])
            matrix[i][j] = max([
                (matrix[i - 1][j - 1][0] + replace, DIAGONAL),
                (matrix[i][j - 1][0] + insert, LEFT),
                (matrix[i - 1][j][0] + delete, UP)
            ])

    # find max score at end of read
    (best_i, best_j, best_score) = 0,0,0
    for i in range(1, len_a + 1):
        #for j in range(1, len_b + 1):
        j = len_b
        if matrix[i][j][0] > best_score:
            best_score = matrix[i][j][0]
            best_i, best_j = i,j
    #print("besti,j",best_i,best_j)
    i, j = best_i, best_j
    align_a = []
    align_b = []

    nb_matches = 0
    nb_columns = 0
    # don't stop until reached beginning of query (or ref)
    while i > 0 and j > 0:
        nb_columns += 1
        if matrix[i][j][1] == DIAGONAL:
            align_a += [a[i - 1]]
            align_b += [b[j - 1]]
            if a[i - 1] == b[j - 1]: nb_matches += 1 
            i -= 1
            j -= 1
        elif matrix[i][j][1] == LEFT:
            align_a += ["-"]
            align_b += [b[j - 1]]
            j -= 1
        else: # UP
            align_a += [a[i - 1]]
            align_b += ["-"]
            i -= 1
    
    blast_identity = 100.0 * nb_matches / nb_columns if nb_columns > 0 else 0
    return best_score, align_a[::-1], align_b[::-1], blast_identity

"""
read1 = reads[0]
print("testing",read1)
res = semiglobal_align(reference,read1)
res = semiglobal_align(reference,read1[::-1])
print(res[0],res[1],res[2])
"""


from multiprocessing import Pool
from contextlib import closing # python/pypy 2 compat,  https://stackoverflow.com/questions/25968518/python-multiprocessing-lib-error-attributeerror-exit
import time

def align(arg):
    read_id, read = arg
    global reference
    fwd = semiglobal_align(reference,read)
    rev = semiglobal_align(reference,read[::-1])
    aln = rev if rev[0] > fwd[0] else fwd
    #time.sleep(0.0001)
    #aln=(0,[0],[0],0)
    aln_score = aln[0]
    # my old identity definiton was:
    #identity = ((100.0*aln_score)/(1.0*len(read)))
    # this is more of a score than an identity.
    # Will now compute BLAST identity, as per
    # https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
    # i.e. number of matches divided by number of columns
    identity = aln[3]
    #print(read,aln_score)
    #print("read identity: %.2f%%" % identity)
    return identity, aln, read_id, read

def process_reads(reads,filename):
    id_dict = dict() # stores identities
    aln_dict = dict() # stores alignments
    orig_dict = dict() # stores original read 'sequences' (of minimizers)

    with closing(Pool(nb_processes)) as p:
        aln_results = p.map(align,reads)

    for (identity, aln, read_id, read) in aln_results:
        id_dict[read_id] = identity
        aln_dict[read_id] = (aln[1],aln[2])
        orig_dict[read_id] = read

    identities = id_dict.values()
    mean_identity = sum(identities) / (1.0*len(identities))
    print("for",filename,"mean read identity: %.2f%%" % mean_identity)
    return id_dict, aln_dict, orig_dict

id_r1, aln_r1, orig_r1 = process_reads(reads, file2)
if reads2 is not None:
    id_r2, aln_r2, orig_r2 = process_reads(reads2, file3)

def short_name(read_id):
    return read_id[:12]+".." if len(read_id) > 12 else read_id

nb_better = 0
nb_nochange = 0
nb_worse  = 0
for seq_id in id_r1:
    if seq_id in id_r2:
        ir1 = id_r1[seq_id]
        ir2 = id_r2[seq_id]
        print("read",short_name(seq_id),"uncor: %0.2f" % ir1,"cor: %0.2f" % ir2)
        if ir1 < ir2:
            nb_better += 1
        elif ir2 < ir1:
            nb_worse += 1
        else:
            nb_nochange += 1
    
        # poa stats
        if file_poa is not None:
            tp, fp, fn = evaluate_poa.eval_poa(seq_id, poa_d, poa_d_itv)
            print("POA retrieval TP: %d  FP: %d  FN: %d" % (tp,fp,fn))

        # jaccard stats
        poa_template = set(orig_r1[seq_id])
        mean_jac = 0
        for poa_seq_id in set(poa_reads[seq_id]):
            poa_r1 = set(orig_r1[poa_seq_id])
            jac_r1 = len(poa_template & poa_r1) / len(poa_template | poa_r1)
            if poa_seq_id in orig_r2:
                poa_r2 = set(orig_r2[poa_seq_id])
                jac_r2 = len(poa_template & poa_r2) / len(poa_template | poa_r2)
            else:
                jac_r2 = -1
            #print("Jac uncor: %.2f    Jac cor: %.2f" % (jac_r1, jac_r2))
            mean_jac += jac_r1
        mean_jac /= len(set(poa_reads[seq_id]))
        print("Mean Jaccard between template and POA reads: %.2f" % mean_jac)

        debug_aln = True 
        if debug_aln:
            print("alignment of uncorrected read",short_name(seq_id)," (len %d) to ref:" % len(orig_r1[seq_id]))
            #print(orig_r1[seq_id]) # print original read sequence of minimizers
            aln = aln_r1[seq_id]
            print("\t".join(map(str,aln[0])))
            print("\t".join(map(str,aln[1])))
            print("and now the corrected read (len %d) alignment:" % len(orig_r2[seq_id]))
            #print(orig_r2[seq_id])
            aln = aln_r2[seq_id]
            print("\t".join(map(str,aln[0])))
            print("\t".join(map(str,aln[1])))

            print("---")

print(nb_better,"reads improved")
print(nb_nochange,"reads unchanged")
print(nb_worse,"reads made worse")
id_r1, aln_r1, orig_r1 = process_reads(reads, file2)
if reads2 is not None:
    id_r2, aln_r2, orig_r2 = process_reads(reads2, file3)
