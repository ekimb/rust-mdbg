def parse(filename):
    # read [.sequences] file
    kmer_to_seq = dict()
    kmer_abundance = dict()
    node_minims = dict()
    k, l = 0, 0
    kmer_origins = dict()
    minim_shift = dict()
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
        minims = tuple(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-5]))
        if spl[-4] != "PLACEHOLDER": 
            abundance = int(spl[-4])
            kmer_abundance[minims] = abundance
        origin = spl[-2]
        seq = spl[-5]
        shifts = tuple(map(lambda x: int(x.strip('(').strip(')').replace(',','')),spl[-2:]))
        node_minims[seq_id] = minims
        kmer_to_seq[minims] = seq
        kmer_origins[minims] = origin
        minim_shift[seq_id] = shifts
    return k, l, node_minims, kmer_to_seq, kmer_abundance, kmer_origins, minim_shift
