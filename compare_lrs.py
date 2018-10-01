from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from random import sample
import logging
from logging import info
from itertools import combinations

logging.basicConfig(filename='info.log',level=logging.DEBUG)

parser = ArgumentParser()

parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("--blacklistfile", help="Blacklist File")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=20)
parser.add_argument("--mincontigs", help="Minimal number of contigs for long read to be considered", type=int, default=2)

args = parser.parse_args()

reads = {}
greads = {}
cgreads = []
contigs = {}


for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

sr_distances = {}
with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        ori1 = sline[0].split("_")[0][0]
        ori2 = sline[0].split("_")[1][0]
        if ori1 != ori2:
            continue
        if ori1 == "-":
            continue
        [ctg1, ctg2] = sline[0].replace("+","").replace("-","").split("_")
        if sline[1] == "NA":
            continue
        if float(sline[2]) < args.mindepth:
            continue
        moddist = float(sline[1])
        # sanity check on the distance
        if moddist + contigs[ctg2] < 0:
            continue
        if ctg1 in sr_distances:
            sr_distances[ctg1][ctg2] = moddist
        else:
            sr_distances[ctg1] = {ctg2: moddist}




#blacklist
blacklist = {}
complete_read = set()
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            idx, ctg =  line.strip().split()[0:2]
            
            if ctg == "all":
                complete_read.add(idx)
            else:
                blacklist[idx] = ctg

print(blacklist)

# nanopore reads
lreads = {}
with open(args.efile) as f:
    for line in f:
        #sline = line.split()
        [rid, ctg, t2, t3, t4, scr, ecr, lenr, strand, scc, ecc, lenc, t12, t13, t14, t15, t16] = line.split()
        data = {"contig":ctg,"strand":int(strand),"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
        if args.blacklistfile:
            if rid in blacklist:
                if blacklist[rid] == ctg:
                    continue
            elif rid in complete_read:
                continue
        if rid in lreads:
            lreads[rid]["maps"].append(data)
        else:
            lreads[rid] = {}
            lreads[rid]["length"] = int(lenr)
            lreads[rid]["maps"] = [data]

# these reads are problematic
def has_contigs_double(lr1):
    seen_ctgs = set()
    for ctg in lr1["maps"]:
        if ctg["contig"] in seen_ctgs:
            return True
        else:
            seen_ctgs.add(ctg["contig"])
    return False

#filter for interesting np reads
greads = {}
for rid,lr in lreads.items():
    if has_contigs_double(lr):
        info(rid + " has contigs double")
        continue
    counter = 0
    for item in lr["maps"]:
        if item["contig"].endswith(args.linename):
            counter +=1
            if counter >= args.mincontigs:
                greads[rid] = lr
                break

lrids = []
# turn reads around if necessary
for rid, lr in greads.items():
    bw = 0
    fw = 0
    lrids.append(rid)
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.linename): 
            if mapping["strand"] == 1:
                bw += 1
            elif mapping["strand"] == 0:
                fw += 1
            else:
                raise ValueError("strand: " + str(mapping["strand"]))
    if bw > fw:
        for mapping in lr["maps"]:
            if mapping["contig"].endswith(args.linename): 
                mapping["strand"] = 1 if mapping["strand"] == 0 else 0
                tmp = mapping["scr"]
                mapping["scr"] = lr["length"] - mapping["ecr"]
                mapping["ecr"] = lr["length"] - tmp
                tmp = mapping["scc"]
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - tmp
        
    # turn around and redefine wrong contigs
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.linename): 
            if mapping["strand"] == 1: #define a new contigname and turn it around
                mapping["contig"] = mapping["contig"] + "rc"
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - mapping["scc"]
                mapping["strand"] = 0

def compare_longreads(lr1, lr2):
    l1c = []
    l2c = []
    for m in lr1["maps"]:
        cn = m["contig"]
        if cn.endswith(args.linename):
            l1c.append(cn)
    for m in lr2["maps"]:
        cn = m["contig"]
        if cn.endswith(args.linename):
            l2c.append(cn)
    common_ctgs = set(l1c).intersection(set(l2c))
    return common_ctgs


def get_contig_info(lr,ctg):
    for maps in lr["maps"]:    
        if maps["contig"] == ctg:
            return maps

def get_distances(lr1,lr2, common_ctgs):
    dists = []
    for ctg in common_ctgs:
        m1 = get_contig_info(lr1,ctg)
        m2 = get_contig_info(lr2,ctg)
        dists.append((m1["scr"]-m1["scc"]) - (m2["scr"]-m2["scc"]))
    return dists

def show_distances(lr1,lr2,lr1id,lr2id, common_ctgs):
    for ctg in common_ctgs:
        m1 = get_contig_info(lr1,ctg)
        m2 = get_contig_info(lr2,ctg)
        dist = ((m1["scr"]-m1["scc"]) - (m2["scr"]-m2["scc"]))
        print(lr1id + " + " + lr2id + " - " + ctg + ": " + str(dist))
        
lr_dists = {}
#initialize matrix
for lid in lrids:
    lr_dists[lid] = {lid:([0],[0])}

# get pair of overlapping read
#lread1, lread2 = sample(list(greads.values()), 2)
#while len(compare_longreads(lread1,lread2)) == 0:
#    lread1, lread2 = sample(list(greads.values()), 2)
#common_ctgs = compare_longreads(lread1,lread2)
all_dists = []

combs = combinations(greads.keys(),2)
for lrs in combs:
    lr1 = greads[lrs[0]]
    lr2 = greads[lrs[1]]
    common_ctgs = compare_longreads(lr1,lr2)
    if len(common_ctgs) > 0:
        dists = get_distances(lr1,lr2,common_ctgs)
        lr_dists[lrs[0]][lrs[1]]=(dists,[])
        stdev = (np.std(dists))
        if stdev > 500:
            show_distances(lr1,lr2,lrs[0],lrs[1],common_ctgs)
            print("-" * 200)
            #print(lrs[0]  + " + " + lrs[1])
        all_dists.append(dists)
        #print(dists)

    # short reads !!
    for m1 in lr1["maps"]:
        ctg1 = m1["contig"]
        for m2 in lr2["maps"]:
            ctg2 = m2["contig"]
            if ctg1 in sr_distances:
                if ctg2 in sr_distances[ctg1]:
                    
                    sr_dist = m1["ecr"] - (contigs[ctg1] - m1["ecc"]) + sr_distances[ctg1][ctg2] - (m2["scr"] - m2["scc"])
                    if lrs[1] in lr_dists[lrs[0]]:
                        lr_dists[lrs[0]][lrs[1]][1].append(sr_dist)
                    else:
                        print(lrs[0] + " + " + lrs[1] + ": " + ctg1 + " " + ctg2)
                        lr_dists[lrs[0]][lrs[1]]= ([],[sr_dist])
                        
    
for lrid,lrdists in lr_dists.items():
    for lrid2,dists in lrdists.items():
        if dists[0] == [] and len(dists[1]) > 1:
            print("interesting: " + str(lrid) + " " + str(lrid2) + " " + str(dists[1]))



#for lrid,dists in lr_dists.items():
#    print("-"*200)
#    print(lrid)
#    print(dists)

        
devs = []
for dists in all_dists:
    stdev = (np.std(dists))
    if stdev > 500:
        print(dists)
    devs.append(stdev)

#print(variances)

# get pairs of overlapping reads
plt.hist(devs,150)
plt.yscale('log', nonposy='clip')
plt.show()
