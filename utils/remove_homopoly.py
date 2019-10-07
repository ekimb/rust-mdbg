# aggressive removal of homopolymers in a FASTA file
import sys
prev_char= ""
for line in open(sys.argv[1]):
    res = ""
    for c in line.strip():
        if c == prev_char: 
            continue
        res += c
        prev_char = c
    print(res)
