import sys
all_genomes = set()
for line in open(sys.argv[1]):
    ls = line.split()
    unitig = ls[0]
    genomes = eval(" ".join(ls[1:]))
    for genome in genomes:
        all_genomes.add(genome.split('.')[0])

for genome in all_genomes:
    print(genome)
