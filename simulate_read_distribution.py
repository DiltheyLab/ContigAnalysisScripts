import numpy as np
from numpy.random import exponential as expo
from random import randint

lchr6 =170805979
smhc =  28510120
emhc =  33480577
lmhc = emhc - smhc

# LSK108
# Number of reads in chr6 and or in contigs 
nreads = 37844 + 1085
nreads = nreads*3
#nreads = int(nreads/2)
good = 0
bad = 0
# the length can be modelled as an exponential distribution
for run in range(0,1000):
    mhc = np.zeros(lmhc)
    draws = expo(8000,[nreads,1])+1200
    print("  " + str(run+1)+"\r" , end='', flush = True)

    for draw in draws:
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
            mhc[bidx:eidx] += 1

    if 0 not in mhc:
        good += 1
    else:
        bad += 1

print("Good: " + str(good) + "\tBad: " + str(bad))
            
