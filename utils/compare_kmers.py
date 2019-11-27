import sys
if len(sys.argv) < 3 or ".sequences" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: [file1.sequences] [file2.sequences]\n will compare the set of kmers in those two files\n")

# tip: make file1 be the .sequences file constructed from a reference genome
# and file2 be one made from the reads
file1= sys.argv[1]
file2= sys.argv[2]

def parse_file(filename):
    res = set()
    for line in open(filename):
        if line.startswith('#'):
            continue 
        spl = line.split()
        unitig_id = spl[0]
        minims = tuple(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-2]))
        seq = spl[-2]
        res.add(minims)
    return res

kmers1 = parse_file(file1)
kmers2 = parse_file(file2)

print(len(kmers1),"kmers in",file1)
print(len(kmers2),"kmers in",file2)

#print("intersection:",len(kmers1 & kmers2))

kmers1_in_kmers2 = len([x for x in kmers1 if x in kmers2])
kmers1_not_in_kmers2  = len(kmers1) - kmers1_in_kmers2
print("number of kmers from %s that are in %s" % (file1,file2), kmers1_in_kmers2, "(%.2f)%%" % ((1.0*kmers1_in_kmers2)/len(kmers1)*100.0),",", kmers1_not_in_kmers2,"are not")
