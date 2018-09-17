import numpy as np
from numpy.random import exponential as expo
from random import randint,sample
from Bio import SeqIO
from argparse import ArgumentParser
from bisect import bisect_left, bisect_right

parser = ArgumentParser() 
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("--lenfile", help="Length of Reads File")
parser.add_argument("--ul_lenfile", help="Length of Ultralong-Reads File")
args = parser.parse_args()

contigs = {}

for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

lengths = []
if args.lenfile:
    with open(args.lenfile) as f:
        for line in f:
            lengths.append(int(line.rstrip()))
ul_lengths = []
if args.ul_lenfile:
    with open(args.ul_lenfile) as f:
        for line in f:
            ul_lengths.append(int(line.rstrip()))


lchr6 =170805979
smhc =  28510120
emhc =  33480577
lmhc = emhc - smhc

pos = 0
nr = 0
crucial_points = []
while pos < lmhc:
    ctgs = sample(contigs.keys(),1)
    crucial_points.append(pos + contigs[ctgs[0]])
    pos += contigs[ctgs[0]]
    #print(pos)
    nr += 1

cp_fixed = {}

# LSK108
# Number of reads in chr6 and or in contigs 
nreads = 37844 + 1085
#nreads = nreads*5
nul_reads = 10511 - 9 + 4735 - 2208
nreads = int(nreads*2)
#nul_reads = nul_reads*2

def find_ge(a, x):
    'Find leftmost index with item greater than or equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        return i
    else:
        return len(a)-1
    raise ValueError
#nreads = int(nreads/2)
good = 0
bad = 0
unbridged = 0
runs = 1000

for run in range(0,runs):
    for cp in crucial_points:
        cp_fixed[cp] = 0
    #mhc = np.zeros(lmhc)
    # the length can be modelled as an exponential distribution
    if args.lenfile:
        draws = sample(lengths,nreads)
    else:
        draws = expo(8000,[nreads,1])+1000
    if args.ul_lenfile:
        ul_draws = sample(ul_lengths,nul_reads)
    else:
        ul_draws = []
    print("  " + str(run+1)+"\r" , end='', flush = True)

    all_draws = draws + ul_draws

    for draw in all_draws:
        beg = randint(1,int(lchr6-draw))
        end = beg+draw
        bidx = 0
        eidx = 0
        if (beg > smhc and beg < emhc) or (end > smhc and end < emhc):
            if beg > smhc:
                bidx = int(beg-smhc -1)
            else:
                bidx = 0
            if end < emhc:
                eidx = int(bidx + draw)
            else:
                eidx = int(lmhc)
            if eidx - bidx < 100:
                continue
            cpil = bisect_right(crucial_points, bidx+50)
            cpir = bisect_right(crucial_points, eidx-50)
            if cpil < cpir:
                for cp in crucial_points[cpil:cpir+1]:
                    cp_fixed[cp] += 1

    #print(cp_fixed)
    

    if 0 not in cp_fixed.values():
        good += 1
    else:
        bad += 1
    unbridged += list(cp_fixed.values()).count(0)

print("Good: " + str(good) + "\tBad: " + str(bad) +"\tUnbridged: " + str(unbridged/runs))
            
