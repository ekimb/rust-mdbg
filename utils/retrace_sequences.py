import sys
if len(sys.argv) < 3 or ".gfa" not in sys.argv[2] or ".sequences" not in sys.argv[1]:
    exit("input: [origin.gfa] [origin.sequences] [target.gfa]\n will propagate info from [origin.sequences] to [target].sequences")

# read [origin.sequences] file
d_minims = dict()
for line in open(sys.argv[1]):
    spl = line.split()
    id = spl[0]
    minims = list(map(lambda x: int(x.strip('[').strip(']').replace(',','')),spl[1:-1]))
    d_minims[id] = minims

def chain_minimizers(info):
    chain = []
    k = len(d_minims[next(iter(d_minims))])
    for (pos, id, ori) in info:
        ms = d_minims[id]
        if len(chain) > 0:
            if chain[-(k-1):] == ms[:k-1]:
                pass
            elif chain[-(k-1):] == ms[::-1][:k-1]:
                ms = ms[::-1]
            else:
                chain = chain[::-1]
                if chain[-(k-1):] == ms[:k-1]:
                    pass
                elif chain[-(k-1):] == ms[::-1][:k-1]:
                    ms = ms[::-1]
                else:
                    print(len(chain),"chain:", chain[-k:], "ms:", ms)
                    print(chain[-(k-1):], ms[::-1][:k-1])
                    exit("unexpected element to chain")
        if len(chain) > 0:
            assert(chain[-(k-1):] == ms[:k-1])
            chain += [ms[-1]]
        else:
            chain = ms
    return chain

output_filename = '.'.join(sys.argv[2].split('.')[:-1])+".sequences"
output = open(output_filename,'w')
def process_unitig(name, info):
    #print("new chain",name,len(info))
    minims = chain_minimizers(info)
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
