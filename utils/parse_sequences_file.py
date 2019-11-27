def parse(filename):
    # read [.sequences] file
    kmer_to_seq = dict()
    kmer_abundance = dict()
    k, l = 0, 0
    for line in open(filename):
        # ignore #'s lines except for getting the k value
        if line.startswith('#'):
            if line.startswith('# k = '):
                k = int(line.split()[-1])
            if line.startswith('# l = '):
                l = int(line.split()[-1])
            continue
        spl = line.split()
        seq_id = spl[0]
        minims = tuple(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-2]))
        abundance = int(spl[-1])
        seq = spl[-2]
        kmer_to_seq[minims] = seq
        kmer_abundance[minims] = abundance
    return k, l, kmer_to_seq, kmer_abundance
