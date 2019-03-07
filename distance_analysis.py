from argparse import ArgumentParser
from operator import itemgetter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from itertools import combinations
from scaffold import Scaffolds


parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("summaryfile", help="Contig distance summary file")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("--include_ambigious", help="Include ambigious contigs", action="store_true", default=False)

args = parser.parse_args()



reads = {}
greads = {}
cgreads = []
ambigious_contigs = set()

blacklist = {}
blacklist_fullread = set()
blacklist_contigs = set()
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist_contigs.add(sline[1])
            if sline[1] == "all":
                blacklist_fullread.add(sline[0])
            else:
                blacklist[sline[0]] = sline[1]

lrs = Scaffolds(args.inputfiles, blacklist, args.linename)
lrs.filter_contigcounts(2)
lrs.turn_longreads_around()
lrs.sort_contigs_in_reads()
lrs = lrs.lreads
            

#print(greads)

distances = {}
# get distances of all neighbouring overlaps
'''
for rid in greads:
    oviter = iter(greads[rid]["overlaps"])
    #print(ov)
    try:
        ovold = next(oviter)
        if ovold["contig"].startswith("chr"):
            continue
    except StopIteration:
        continue
    while True:
        try:
            ovnew = next(oviter)
            if ovnew["contig"].startswith("chr"):
                continue
        except StopIteration:
            break
        #print("distance between " + ovold["contig"] + " and " + ovnew["contig"] + ": " + str(ovnew["scr"] - ovold["ecr"]- ovold["lc"] + ovold["ecc"] - ovnew["scc"]))
        if ovnew["contig"] == ovold["contig"]:
            continue
        if ovnew["strand"] == ovold["strand"]:
            distance = ovnew["scr"] - ovold["ecr"]- (ovold["lc"] - ovold["ecc"]) - ovnew["scc"] + 1
        else:
            continue
        if int(ovold["contig"].rstrip(args.linename)) < int(ovnew["contig"].rstrip(args.linename)):
            cstring = ovold["contig"] + "_" + ovnew["contig"]
        else:
            cstring = ovnew["contig"] + "_" + ovold["contig"]
    
        if cstring in distances:
            distances[cstring].append(distance)
        else:
            distances[cstring] = [distance]
        ovold = ovnew
'''

# get distances of all overlaps
for rid in lrs:
    for item in combinations(lrs[rid]['maps'], 2):
        ovold, ovnew = item
        if ovnew["name"].startswith("chr") or ovold["name"].startswith("chr"):
            continue
        if ovnew["name"] in ambigious_contigs or ovold["name"] in ambigious_contigs:
            continue
        if "_" in ovnew["name"] or "_" in ovold["name"]:
            continue
        if ovnew["name"] == ovold["name"]:
            continue
        if ovnew["strand"] == 1 or ovold["strand"] == 1:
            continue
        distance = ovnew["scr"] - ovold["ecr"]- (ovold["lenc"] - ovold["ecc"]) - ovnew["scc"] + 1
        cstring = ovold["name"] + "_" + ovnew["name"]
        if cstring in distances:
            distances[cstring].append(distance)
        else:
            distances[cstring] = [distance]
        #print(combo)
    

        
for key, value in distances.items():
    #print(str(key) + " " + str(value)) 
    if len(value)>1:
        #print(str(key) + "\t" + str(value)) 
        pass

distances2 = {}

with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        ctg1 = sline[0].split("_")[0].strip("+").strip("-")
        ctg2 = sline[0].split("_")[1].strip("+").strip("-")
        if line.startswith("-"):
            continue
        if sline[1] == "NA":
            continue
        if float(sline[4]) > 2:
            continue
        moddist = float(sline[1])
        #if int(ctg1.rstrip(args.linename)) < int(ctg2.rstrip(args.linename)):
        cstr = ctg1+"_"+ctg2
        #else:
        #    cstr = ctg2+"_"+ctg1
        if cstr in distances2:
            if abs(moddist) < abs(distances2[cstr]):
                distances2[cstr] = moddist
        else:    
            distances2[cstr] = moddist

for name, dist in distances.items():
    if name in distances2:
        dist2 = distances2[name]
    else:
        dist2 = "-"
    name1, name2 = name.split("_")
    print("\t".join([name1, name2, str(dist), str(dist2)]))
            
        
df = pd.DataFrame.from_dict([distances, distances2])
#df.rename(index=
dc = df.T.rename(columns={0:'longread',1:'shortread'})
dc["longread_mean"] = dc.longread.apply(np.mean)

#dc['longread']  = np.mean(dc.longread)
dd = dc.dropna()

#get interesting differences
#print(dd[abs(dd['longread_mean'] - dd['shortread']) > 150])
#print(dd)
sthsth = []
for item in dd['longread_mean']:
    sthsth.append(item < 0)
for idx, item in enumerate(dd['shortread']):
    sthsth[idx] = sthsth[idx] or item < 0
#for name in dd[sthsth].index.values:
#    print(name)
#print(dd[dd['longread_mean'] <= -20])
#print(dd.index.values)


plt.scatter(dd['longread_mean'], dd['shortread'],s= 6, alpha = 0.7)
plt.xlabel("Long Read Distances (mean: " + "{:.3f}".format(np.mean(dd['longread_mean'])) + ")")
#plt.xlabel("Long Read Distances")
plt.ylabel("Short Read Distances (mean: " + "{:.3f}".format(np.mean(dd['shortread'])) + ")")
#plt.ylabel("Short Read Distances")
plt.savefig('distances_scatter.pdf')
