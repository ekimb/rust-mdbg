import sys
import glob
import re

# usage: find -name "*.n50"|xargs cat| python make_table.py
#
#or, when dealing with subfolders having the param names, and the .n50 file is meaningless:
# find -name "*.n50"|xargs tail -n +1| python make_table.py

#folder="/pasteur/scratch/public/rchikhi/mouse-mdbg/"

print("cvg,k,l,d,n50")


"""
assembly:/pasteur/scratch/public/rchikhi/mouse-mdbg/mouse-K11-L10-D0.001.msimpl.fa
number of contigs/scaffolds:10883
assembly size:245671386
largest contig/scaffold:596186
N50:219897
N90:10960
"""

cvg=None
for line in sys.stdin:
            line = line.strip()
            line = line.strip('.n50')
            if 'assembly:' in line:
                numbers = re.findall('[0-9]+', line)
                #print(line,numbers)
                if len(numbers) == 4:
                    k,l,osef,d=numbers
                    cvg=""
                elif len(numbers) == 6:
                    osef,cvg,k,l,osef,d=numbers
                    cvg=float("0."+cvg)
                elif len(numbers) == 5: # when coverage=1
                    cvg,k,l,osef,d=numbers
                    assert(float(cvg)==1)
                    cvg=1
                else:
                    k=numbers[0]
                k=int(k)
                try:
                    d=float("0."+d)
                except:
                    print("error processing line: " + line + " with recorded params " + str([cvg,k,l,d]))
                l=int(l)
                if l not in [10,11,12,13,14,15]:
                    exit("bad format")
                #print(line,k,l,d)
            if "N50" in line:
                n50 = line.split(':')[-1]
                print(cvg,k,l,d,n50,sep=",")
                numbers = re.findall('[0-9]+', line)
            if "==>" in line:
                numbers = re.findall('[0-9]+', line)
                if len(numbers) == 5:
                    osef,d,l,k,osef = numbers

