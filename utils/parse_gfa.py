def parse(filename):
    # read [.gfa] file, just get abundances for now
    kmer_abundance = dict()
    for line in open(filename):
        if line.startswith('S'):
            #S       840     *       LN:i:1  KC:i:1
            spl = line.split()
            seq_id = spl[1]
            for field in spl:
                if field.startswith("KC:i"):
                    abundance = int(field.split(':')[-1])
                    kmer_abundance[seq_id] = abundance
    return kmer_abundance
