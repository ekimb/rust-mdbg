import sys
if len(sys.argv) < 3 or ".gfa" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: [origin.sequences] [target.gfa]\n will propagate info from [origin.sequences] to [target].sequences")

# read [origin.sequences] file
d_minims = dict()
from parse_sequences_file import parse
k, l, node_minims, kmer_seq, kmer_abundance, origins1 = parse(sys.argv[1])
d_minims = node_minims

def chain_minimizers(info, unitig_name): # unitig_name is just for debug
    abunds = []
    chain = []
    for (chain_number,(pos, node_id, ori)) in enumerate(info):
        # FIXME for some reason I didn't use the 'ori' field but it could actually help
        ms = d_minims[node_id]
        abund = kmer_abundance[ms]
        if len(chain) > 0:
            if chain[-(k-1):] == ms[:k-1]:
                pass
            elif chain[-(k-1):] == ms[::-1][:k-1]:
                ms = ms[::-1]
            else:
                bad = False
                if chain_number == 1: # can only reverse the first element of the chain during attempt to chain the the second element
                    chain = chain[::-1]
                    if chain[-(k-1):] == ms[:k-1]:
                        pass
                    elif chain[-(k-1):] == ms[::-1][:k-1]:
                        ms = ms[::-1]
                    else: 
                        bad = True
                else:
                    bad = True
                if bad:
                    # some extensive debugging information
                    print("chain (size %d):" % len(chain), "last k=%d elements:" % k,chain[-k:])
                    print("                 first k=%d elements:" % k,chain[:k])
                    print("to be chained with node id %s (size %d):" %(node_id,len(ms)), ms)
                    if chain == ms or chain == ms[::-1]:
                        print("!!warning!! chain == ms or chain == ms[::-1]")
                    print(", so tested overlap:",chain[-(k-1):])
                    print("        with either:",ms[:k-1])
                    print("                 or:",ms[::-1][:k-1])
                    continue
                    #exit("unexpected element to chain (unitig %s)" % unitig_name)
        if len(chain) > 0:
            assert(chain[-(k-1):] == ms[:k-1])
            chain += ms[k-1:][::]
            print("ID %s abund %d" % (node_id, abund))
            abunds.append(abund)

        else:
            print("ID %s abund %d" % (node_id, abund))
            abunds.append(abund)
            chain = ms[::]
            # small note to myself:
            # gfaview re-uses unitig's across simplifications. i.e. the same origin unitig may be found in two different 'a' lines
    return chain, abunds

output_filename = '.'.join(sys.argv[2].split('.')[:-1])+".sequences"
output = open(output_filename,'w')
output.write("# k = %d\n" % k)
output.write("# l = %d\n" % l)
def process_unitig(name, info):
    MIN_ABUNDANCE = 3
    print("new chain",name,"len",len(info),"contents:",info)
    minims, abunds = chain_minimizers(info, name)
    passed = [x for x in abunds if x > MIN_ABUNDANCE]
    if len(passed) != 0:
        print("Passed nodes %d" % len(passed))
        output.write("%s\t%s\tPLACEHOLDER\tPLACEHOLDER\tPLACEHOLDER\n"% (name,minims))

# read [target.gfa] file
current_unitig_name = ""
current_unitig_info = []
for line in open(sys.argv[2]):
    if not line.startswith('a'): continue
    # a       utg0010623      0       490197  +       100
    spl = line.split()
    unitig_name = spl[1]
    unitig_pos = spl[2]
    seq_id = spl[3]
    ori = spl[4]
    if unitig_name != current_unitig_name:
        if current_unitig_name != "":
            process_unitig(current_unitig_name, current_unitig_info)
        current_unitig_name = unitig_name
        current_unitig_info = []
    current_unitig_info += [(unitig_pos, seq_id, ori)]
process_unitig(current_unitig_name, current_unitig_info)
