# this script artificially breaks self-loop nodes by cutting one of the extremity edges
# e.g. 
# L x + y - 10M
# L x + y + 10M
# will remove one of the two lines

seen_edges = set()
import sys, traceback
for line in open(sys.argv[1]):
    if not line.startswith('L'): 
        print(line.strip())
        continue
    try:
        e_in, e_out = line.split()[1], line.split()[3]
        e_tuple = tuple(sorted([e_in, e_out]))
        to_remove = e_tuple in seen_edges
        seen_edges.add(e_tuple)
        if not to_remove:
            print(line.strip())
    except:
        print("at line:",line.strip())
        traceback.print_exc(file=sys.stdout)
        exit(1)
