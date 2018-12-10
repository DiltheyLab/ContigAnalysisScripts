from argparse import ArgumentParser
from operator import itemgetter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from itertools import combinations


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="Contig distance summary file")
parser.add_argument("cellline", help="Name of cell line")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("--include_ambigious", help="Include ambigious contigs", action="store_true", default=False)

args = parser.parse_args()


reads = {}
greads = {}
cgreads = []

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
ambigious_contigs = set()

# nanopore read info
with open(args.efile) as f:
    for line in f:
        sline = line.split()
        rid = sline[0]
        ctg = sline[1]
        strand = int(sline[8])
        scr = int(sline[5])
        ecr = int(sline[6])
        scc = int(sline[9])
        ecc = int(sline[10])
        lc = int(sline[11])
        payload = {"contig":ctg,"strand":strand,"scr":scr,"ecr":ecr,"scc":scc,"ecc":ecc,"lc":lc}
        if ctg.startswith("b_") or ctg.startswith("c_"):
            ambigious_contigs.add(ctg.split("_")[1])
        if ctg in blacklist_contigs:
            continue
        elif rid in blacklist:
            if blacklist[rid] == ctg:
                continue
        elif rid in blacklist_fullread:
            continue

        if rid in reads:
            reads[rid]["overlaps"].append(payload)
        else:
            reads[rid] = {}
            reads[rid]["length"] = int(sline[7])
            reads[rid]["overlaps"] = [payload]

# get interesting reads
greadst = {}
for rid in reads:
    counter = 0
    for item in reads[rid]["overlaps"]:
        if item["contig"].endswith(args.cellline):
            counter += 1
    if counter >= 2:
        greadst[rid] = reads[rid]
            
# sort contigs by left coordinate
for rid in greadst:
    #print(reads[rid]["overlaps"])
    soverlaps = sorted(reads[rid]["overlaps"], key = itemgetter("scr"))
    greads[rid] = greadst[rid]
    greads[rid]["overlaps"]=soverlaps

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
        if int(ovold["contig"].rstrip(args.cellline)) < int(ovnew["contig"].rstrip(args.cellline)):
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
for rid in greads:
    for combo in combinations(greads[rid]['overlaps'], 2):
        ovnew = combo[1]
        ovold = combo[0]
        if ovnew["contig"].startswith("chr") or ovold["contig"].startswith("chr"):
            continue
        if ovnew["contig"] in ambigious_contigs or ovold["contig"] in ambigious_contigs:
            continue
        if "_" in ovnew["contig"] or "_" in ovold["contig"]:
            continue
        if ovnew["contig"] == ovold["contig"]:
            continue
        if ovnew["strand"] == ovold["strand"]:
            distance = ovnew["scr"] - ovold["ecr"]- (ovold["lc"] - ovold["ecc"]) - ovnew["scc"] + 1
        else:
            continue
        #if abs(distance) > 2000:
        #    continue
        if int(ovold["contig"].rstrip(args.cellline)) < int(ovnew["contig"].rstrip(args.cellline)):
            cstring = ovold["contig"] + "_" + ovnew["contig"]
        else:
            cstring = ovnew["contig"] + "_" + ovold["contig"]
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
        if sline[1] == "NA":
            continue
        if float(sline[4]) > 2:
            continue
        moddist = float(sline[1])
        if int(ctg1.rstrip(args.cellline)) < int(ctg2.rstrip(args.cellline)):
            cstr = ctg1+"_"+ctg2
        else:
            cstr = ctg2+"_"+ctg1
        if cstr in distances2:
            if abs(moddist) < abs(distances2[cstr]):
                distances2[cstr] = moddist
        else:    
            distances2[cstr] = moddist
            
        
df = pd.DataFrame.from_dict([distances, distances2])
#df.rename(index=
dc = df.T.rename(columns={0:'longread',1:'shortread'})
dc["longread_mean"] = dc.longread.apply(np.mean)

#dc['longread']  = np.mean(dc.longread)
dd = dc.dropna()

#get interesting differences
print(dd[abs(dd['longread_mean'] - dd['shortread']) > 150])


plt.scatter(dd['longread_mean'], dd['shortread'],s= 6, alpha = 0.4)
plt.xlabel("Long Read Distances (mean: " + "{:.3f}".format(np.mean(dd['longread_mean'])) + ")")
#plt.xlabel("Long Read Distances")
plt.ylabel("Short Read Distances (mean: " + "{:.3f}".format(np.mean(dd['shortread'])) + ")")
#plt.ylabel("Short Read Distances")
plt.savefig('distances_scatter.pdf')
