from collections import defaultdict, Counter
import sys
if len(sys.argv) < 3 or ".ec_data" not in sys.argv[2] or ".ec_data" not in sys.argv[1]:
    exit("input: [reference.ec_data] [reads.sequences]\n will attempt EC on reads and use reference to tell if it went ok\n")

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

def normalize(t):
    return min(t,t[::-1])

lmers = Counter()
success_dict = defaultdict(Counter)
l=3
for read in reads:
    for i in range(len(read)-l+1):
        lmer = normalize(tuple(read[i:i+l]))
        lmers[lmer] += 1
        if i < len(read)-l:
            success_dict[lmer[:-1]][lmer[-1]] += 1
print(lmers)

for read in reads:
    print(read)
    for i in range(len(read)-l+1):
        lmer = normalize(tuple(read[i:i+l]))
        if lmers[lmer] < 20:
            print(i,"weak, successors",success_dict[lmer[1:]])
