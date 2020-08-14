"""
 this script takes as input:
 - the .poa.ec_data file for a set of reads, in the format: template \t read1 \t read2 \t etc..
"""

min_overlap = 1000


from intervaltree import Interval, IntervalTree

def parse_file(filename):
    d = dict()
    poa_reads = dict()
    for line in open(filename):
        template = line.split()[0]
        d[template] = []
        poa_reads[template] = []
        for read in line.split()[1:]:
            #SYN_658_52663_60655_0_-_60663_1_._NC{004353.4$Drosophila$melanogaster$chromosome$4
            read_start = int(read.split('_')[2])
            read_end   = int(read.split('_')[3])
            d[template] += [(read_start,read_end)]
            poa_reads[template] += [read]    
        if (len(d[template]) != len(set(d[template]))):
            print("warning, duplicated reads retrieved for POA template",template,":",len(d[template]),"!=",len(set(d[template])))
    return d, poa_reads

# just a function that computes the overlap between two intervals.. hopefully it's not buggy
def overlap(itv, s, e):
    a,b = itv.begin, itv.end
    if a <= s:
        # a---b
        #   s---e or
        # a-------b
        return min(b-s,e-s) if b > s else 0
    else:
        #  a---b
        # s---e or
        #  a-----b
        return min(e-a,b-a) if e > a else 0

def get_intervals(start,end,it):
    all_intervals = it[start:end]
    # do a filter for long enough overlaps
    return set([(itv.begin,itv.end) for itv in all_intervals if overlap(itv,start,end) > min_overlap])

def prepare_eval_poa(filename):
    d, poa_reads = parse_file(filename)

    whole_it = IntervalTree(set([Interval(*item) for sublist in d.values() for item in sublist]))

    return d, whole_it, poa_reads

def eval_poa(template, d, d_itv):
    template_start = int(template.split('_')[2])
    template_end   = int(template.split('_')[3])

    truth_intervals = get_intervals(template_start, template_end, d_itv)
    poa_intervals   = set(d[template])

    tp = len([x for x in poa_intervals if x in     truth_intervals])
    fp = len([x for x in poa_intervals if x not in truth_intervals])
    fn = len([x for x in truth_intervals if x not in poa_intervals])

    return tp, fp, fn


if __name__ == "__main__": 
    import sys
    if len(sys.argv) < 2 or ".ec_data" not in sys.argv[1]:
        exit("input: [reads.ec_data] [min_overlap=2000]\n will evaluate if POA retrieves the right set of reads, those that overlap with the template over min_overlap base \n")

    filename = sys.argv[1]
    if len(sys.argv) > 2:
        min_overlap = int(sys.argv[2])
    
    d, d_itv = prepare_eval_poa(filename)
    
    for template in d:
        tp, fp, fn = eval_poa(template, d, d_itv)
        print("TP: %d  FP: %d  FN: %d" % (tp,fp,fn))


