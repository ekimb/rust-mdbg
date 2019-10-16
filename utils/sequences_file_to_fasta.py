import sys
if len(sys.argv) < 3 or ".fa" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: [final.sequences] [final.fasta]\n will take the sequences column and put it in fasta format, that's all\n")

output = open(sys.argv[2],"w")
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue 
    spl = line.split()
    unitig_id = spl[0]
    minims = tuple(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-2]))
    seq = spl[-2]
    output.write(">%s\n%s\n" % (unitig_id,seq))
output.close()
     
