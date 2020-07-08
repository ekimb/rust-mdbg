import sys
if len(sys.argv) < 3 or ".gfa" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: [origin.sequences] [target.gfa]\n will propagate info from [origin.sequences] to [target].sequences")

# read [origin.sequences] file
d_minims = dict()
k, l = 0, 0
for line in open(sys.argv[1]):
    # ignore #'s lines except for getting the k value
    if line.startswith('#'):
        if line.startswith('# k = '):
            k = int(line.split()[-1])
        if line.startswith('# l = '):
            l = int(line.split()[-1])
        continue
    spl = line.split()
    node_id = spl[0]
    minims = list(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-2]))
    d_minims[node_id] = minims

def chain_minimizers(info, unitig_name): # unitig_name is just for debug
    chain = []
    for (chain_number,(pos, node_id, ori)) in enumerate(info):
        ms = d_minims[node_id]
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
                    print("to be chained with node id %s (size %d):" %(node_id,len(ms)), ms)
                    if chain == ms or chain == ms[::-1]:
                        print("!!warning!! chain == ms or chain == ms[::-1]")
                    print(", so tested overlap:",chain[-(k-1):])
                    print("        with either:",ms[:k-1])
                    print("                 or:",ms[::-1][:k-1])
                    exit("unexpected element to chain (unitig %s)" % unitig_name)
        if len(chain) > 0:
            assert(chain[-(k-1):] == ms[:k-1])
            chain += ms[k-1:][::]
        else:
            chain = ms[::]
            # small note to myself:
            # gfaview re-uses unitig's across simplifications. i.e. the same origin unitig may be found in two different 'a' lines
    return chain

output_filename = '.'.join(sys.argv[2].split('.')[:-1])+".sequences"
output = open(output_filename,'w')
output.write("# k = %d\n" % k)
output.write("# l = %d\n" % l)
def process_unitig(name, info):
    #print("new chain",name,len(info))
    minims = chain_minimizers(info, name)
    output.write("%s\t%s\tPLACEHOLDER\n"% (name,minims))

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
