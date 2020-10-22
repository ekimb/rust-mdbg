import sys
from parse_sequences_file import parse
k, l, node_minims, kmer_seq, kmer_abundance, origins1, minim_shift = parse(sys.argv[1])

def find_overlap(source, sink):
    if source[1] == "+":
        shift = source[3][0]
    else:
        shift = source[3][1]
    assert(shift < len(source[2]))
    return shift


output_filename = '.'.join(sys.argv[1].split('.')[:-1])+".complete.gfa"
output = open(output_filename,'w')
output.write("H\tVN:Z:1\n")
for line in open(sys.argv[2]):
    if not line.startswith('L'): continue
    #L	2377	+	5976	+	0M
    spl = line.split()
    source_minims = node_minims[spl[1]]
    sink_minims = node_minims[spl[3]]
    source = (spl[1], spl[2], kmer_seq[source_minims], minim_shift[spl[1]])
    sink = (spl[3], spl[4], kmer_seq[sink_minims], minim_shift[spl[3]])
    shift = find_overlap(source, sink)
    overlap_length = len(source[2]) - shift
    overlap_length = min(overlap_length, len(sink[2])-1)

    output.write("S\t%s\t%s\tLN:i:%d\tKC:i:%d\n" % (source[0], source[2], len(source[2]), int(kmer_abundance[source_minims])))
    output.write("S\t%s\t%s\tLN:i:%d\tKC:i:%d\n" % (sink[0], sink[2], len(sink[2]), int(kmer_abundance[sink_minims])))
    output.write("L\t%s\t%s\t%s\t%s\t%dM\n"% (source[0], source[1], sink[0], sink[1], int(overlap_length)))
