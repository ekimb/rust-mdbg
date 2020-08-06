"""
 this script takes as input the .ec_data file for a reference genome (processed by rust-mhdbg with min-abundance of 1)
 as well as the .ec_data file for a set of reads (also processed by rust-mhdbg but with the desired error correction method)
 and outputs the needleman-wunch alignment of each read to the reference (semi-global aln)  
"""
import sys
if len(sys.argv) < 3 or ".ec_data" not in sys.argv[2] or ".ec_data" not in sys.argv[1]:
    exit("input: [reference.ec_data] [reads.ec_data]\n will evaluate accuracy of minimizers in reads\n")

# tip: make file1 be the .ec_data file constructed from a reference genome
# and file2 be one made from the reads
file1= sys.argv[1]
file2= sys.argv[2]

def parse_file(filename):
    res = []
    counter = 0
    for line in open(filename):
        if counter % 4 != 1:
            counter += 1 
            continue
        spl = line.split()
        minims = list(map(int,spl))
        res += [minims]
        counter += 1
    return res


reference = parse_file(file1)
assert(len(reference)==1)
reads = parse_file(file2)

print("loaded",len(reference),"reference,",len(reads),"reads")

reference = reference[0]

# adapted from NW code here https://stackoverflow.com/questions/2718809/how-to-diff-align-python-lists-using-arbitrary-matching-function
# 
def semiglobal_align(a, b):
    # my custom scoring
    # not quite edit distance: match would need to be 0 then 
    # but i dunno how to compute semiglobal edit distance so i'll take that for now!
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

    # don't stop until reached beginning of query (or ref)
    while i > 0 and j > 0:
        if matrix[i][j][1] == DIAGONAL:
            align_a += [a[i - 1]]
            align_b += [b[j - 1]]
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

    return best_score, align_a[::-1], align_b[::-1]

"""
read1 = reads[0]
print("testing",read1)
res = semiglobal_align(reference,read1)
res = semiglobal_align(reference,read1[::-1])
print(res[0],res[1],res[2])
"""

identities = []
for read in reads:
    fwd = semiglobal_align(reference,read)
    rev = semiglobal_align(reference,read[::-1])
    aln = rev if rev[0] > fwd[0] else fwd
    aln_score = aln[0]
    identity = ((100.0*aln_score)/(1.0*len(read)))
    identities += [identity]
    #print(read)
    print("read identity: %.2f%%" % identity)

mean_identity = sum(identities) / (1.0*len(identities))
print("mean read identity: %.2f%%" % mean_identity)
