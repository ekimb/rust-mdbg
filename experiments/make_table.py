import sys
import glob
import re

folder="/pasteur/scratch/public/rchikhi/mouse-mdbg/"
          
print("k l d n50")


if len(sys.argv) >= 2:
    folder=sys.argv[1]
    for n50 in glob.glob(folder+"/*.n50"):
        for line in open(n50):
            """
assembly:/pasteur/scratch/public/rchikhi/mouse-mdbg/mouse-K11-L10-D0.001.msimpl.fa
number of contigs/scaffolds:10883
assembly size:245671386
largest contig/scaffold:596186
N50:219897
N90:10960
"""
            line = line.strip()
            if 'assembly:' in line:
                numbers = re.findall('[0-9]+', line)
                #print(line,numbers)
                k,l,osef,d=map(int,numbers)
                d=float("0."+numbers[-1])
                #print(line,k,l,d)
            if "N50" in line:
                n50 = line.split(':')[-1]
                print(k,l,d,n50)
                
