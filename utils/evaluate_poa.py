"""
 this script takes as input:
 - the .poa.ec_data file for a set of reads, in the format: template \t read1 \t read2 \t etc..
"""

min_overlap = 1000


from intervaltree import Interval, IntervalTree

def parse_file(filename, only_those_reads=None):
    d = dict()
    poa_reads = dict()
    all_reads = dict()
    for line in open(filename):
        template = line.split()[0] 
        template_start = int(template.split('_')[2])
        template_end   = int(template.split('_')[3])
        all_reads[template] = (template_start,template_end) # assumption : all reads are seen as template
        d[template] = []
        if only_those_reads is not None and template not in only_those_reads:
            continue
        poa_reads[template] = []
        for read in line.split()[1:]:
            #SYN_658_52663_60655_0_-_60663_1_._NC{004353.4$Drosophila$melanogaster$chromosome$4
            assert(read.startswith("SYN_"))
            read_start = int(read.split('_')[2])
            read_end   = int(read.split('_')[3])
            d[template] += [Interval(read_start,read_end,read)]
            poa_reads[template] += [read]    
            if read not in all_reads:
                all_reads[read] = (read_start,read_end)
        if (len(d[template]) != len(set(d[template]))):
            print("warning, duplicated reads retrieved for POA template",template,":",len(d[template]),"!=",len(set(d[template])))
    print("parsed",filename,"found",len(all_reads),"reads as templates")
    return d, poa_reads, all_reads

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
    return set([itv for itv in all_intervals if overlap(itv,start,end) > min_overlap])

def prepare_eval_poa(filename,only_those_reads=None):
    d, poa_reads, all_reads = parse_file(filename,only_those_reads)
    
    whole_it = IntervalTree([Interval(*all_reads[read_id], read_id) for read_id in all_reads])

    return d, whole_it, poa_reads

def eval_poa(template, d, d_itv):
    if template not in d:
        print("read",template,"isn't a POA template")
        return [],[],[]

    template_start = int(template.split('_')[2])
    template_end   = int(template.split('_')[3])

    truth_intervals = get_intervals(template_start, template_end, d_itv)
    poa_intervals   = set(d[template])

    tp = [x.data for x in poa_intervals if x in     truth_intervals]
    fp = [x.data for x in poa_intervals if x not in truth_intervals]
    fn = [x.data for x in truth_intervals if x not in poa_intervals]

    return tp, fp, fn


if __name__ == "__main__": 
    import sys
    if len(sys.argv) < 2 or ".ec_data" not in sys.argv[1]:
        exit("input: [reads.ec_data] [min_overlap=2000]\n will evaluate if POA retrieves the right set of reads, those that overlap with the template over min_overlap base \n")

    filename = sys.argv[1]
    if len(sys.argv) > 2:
        min_overlap = int(sys.argv[2])
    
    d, d_itv, poa_reads = prepare_eval_poa(filename)
    
    for template in d:
        tp, fp, fn = eval_poa(template, d, d_itv)
        print("TP: %d  FP: %d  FN: %d" % (len(tp),len(fp),len(fn)))


