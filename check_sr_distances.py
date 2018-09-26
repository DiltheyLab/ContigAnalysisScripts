from argparse import ArgumentParser
from operator import itemgetter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from itertools import combinations


parser = ArgumentParser()
parser.add_argument("summaryfile", help="Contig distance summary file")
parser.add_argument("cellline", help="Name of cell line")
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=25)

args = parser.parse_args()


distances2 = {}

with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        ctg1 = sline[0].split("_")[0]
        ctg2 = sline[0].split("_")[1]
        if sline[1] == "NA":
            continue
        if float(sline[2]) < args.mindepth:
            continue
        moddist = float(sline[1])
        depth = int(sline[2])
        if ctg1.rstrip(args.cellline) < ctg2.rstrip(args.cellline):
            cstr = ctg1+"_"+ctg2
        else:
            cstr = ctg2+"_"+ctg1
        if cstr in distances2:
            if abs(moddist) < abs(distances2[cstr][0]):
                distances2[cstr] = (moddist,depth)
        else:    
            distances2[cstr] = (moddist,depth)
            
dists = []
depths = []


for key,[dist,depth] in distances2.items():
    if dist > -500:
        dists.append(dist)
        depths.append(depth)

print(len(dists))
print(len(depths))
a,bin_edges = np.histogram(dists,bins=range(-500,610,10),weights=depths)
b,bin_edges = np.histogram(dists,bins=range(-500,610,10))

averages=[]

for idx, item in enumerate(a):
    averages.append(a[idx]/b[idx])
    
print(averages)
#plt.subplot(2,2,1)
#plt.hist(dists,bins=range(-500,610,10),weights=depths)
#plt.subplot(2,2,2)
#plt.hist(dists,bins=range(-500,610,10))
#plt.subplot(2,2,3)
plt.plot(range(-500,600,10),averages)
plt.show()
