from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from random import sample
import logging
from logging import info
from itertools import combinations
from collections import defaultdict
import os.path

logging.basicConfig(filename='info.log',level=logging.DEBUG)

parser = ArgumentParser()

parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("--blacklistfile", help="Blacklist File")
parser.add_argument("--overwrite", help="Overwrite preexisting distance matrix.", action = "store_true")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=20)
parser.add_argument("--mincontigs", help="Minimal number of contigs for long read to be considered", type=int, default=2)
parser.add_argument("--unwanted_contigs", help="If given the contigs in this file will not be considered.")

args = parser.parse_args()

reads = {}
greads = {}
cgreads = []
contigs = {}


for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

#unwanted contigs
unwanted_contigs = set()
if args.unwanted_contigs:
    with open(args.unwanted_contigs) as f:
        for line in f:
            ctg =  line.strip()
            unwanted_contigs.add(ctg)
            unwanted_contigs.add("b_" + ctg) # usual prefix for repeated contigs
            unwanted_contigs.add("c_" + ctg)

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
        if ctg1 in unwanted_contigs or ctg2 in unwanted_contigs:
            continue
        moddist = float(sline[1])
        # sanity check on the distance
        if moddist + contigs[ctg1] < 0 or moddist + contigs[ctg2] < 0 :
            continue
        if ctg1 in sr_distances:
            sr_distances[ctg1][ctg2] = moddist
        else:
            sr_distances[ctg1] = {ctg2: moddist}

# some sr distances suck
'''
sr_distances["1046APD"].pop("1550APD")
sr_distances["1488APD"].pop("1553APD")
sr_distances["169APD"].pop("530APD")
sr_distances["2137APD"].pop("530APD")
sr_distances["367APD"].pop("2038APD")
sr_distances["367APD"].pop("398APD")
sr_distances["544APD"].pop("1923APD")
sr_distances["582APD"].pop("635APD")
'''

#blacklist of long reads
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


dm_path = "/home/houwaart/Projects/ImmunoPore/APD/distance_matrix"
# if distance matrix exists there is no need to calculate
if args.overwrite or not os.path.isfile(dm_path):
                
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
            if args.unwanted_contigs:
                if ctg in unwanted_contigs:
                    continue
            if rid in lreads:
                lreads[rid]["maps"].append(data)
                if int(ecr) > lreads[rid]["rm_ecr"]:
                    lreads[rid]["rm_ecr"] = int(ecr)
                if int(scr) < lreads[rid]["lm_scr"]:
                    lreads[rid]["lm_scr"] = int(scr)
            else:
                lreads[rid] = {}
                lreads[rid]["length"] = int(lenr)
                lreads[rid]["maps"] = [data]
                lreads[rid]["rm_ecr"] = int(ecr)
                lreads[rid]["lm_scr"] = int(scr)

    # these reads are problematic
    def has_contigs_double(lr1):
        seen_ctgs = set()
        for ctg in lr1["maps"]:
            if ctg["contig"] in seen_ctgs:
                return ctg
            else:
                if ctg["contig"].endswith(args.linename):
                    seen_ctgs.add(ctg["contig"])
        return ""

    #filter for interesting np reads
    greads = {}
    for rid,lr in lreads.items():
        double_ctg =  has_contigs_double(lr)
        if double_ctg:
            info(rid + " has contigs double " + str(double_ctg))
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
            #print(lr1id + " + " + lr2id + " - " + ctg + ": " + str(dist))
            
    #lr_dist = defaultdict(lambda:([],[]))
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

    for lrs in combinations(greads.keys(), 2):
        lr1 = greads[lrs[0]]
        lr2 = greads[lrs[1]]
        common_ctgs = compare_longreads(lr1, lr2)
        if len(common_ctgs) > 0:
            dists = get_distances(lr1, lr2, common_ctgs)
            lr_dists[lrs[0]][lrs[1]]=(dists, [])
            ndists =[]
            for d in dists:
                ndists.append(-d)
            lr_dists[lrs[1]][lrs[0]] = (ndists, [])
            stdev = (np.std(dists))
            if stdev < 50:
                show_distances(lr1, lr2, lrs[0], lrs[1], common_ctgs)
                #print("-" * 200)
                #print(lrs[0]  + " + " + lrs[1])
            all_dists.append(dists)
            #print(dists)
            

        # short reads !
        for m1 in lr1["maps"]:
            ctg1 = m1["contig"]
            for m2 in lr2["maps"]:
                ctg2 = m2["contig"]
                if ctg1 in sr_distances:
                    if ctg2 in sr_distances[ctg1]:
                        sr_dist = m1["ecr"] + (contigs[ctg1] - m1["ecc"]) + sr_distances[ctg1][ctg2] - (m2["scr"] - m2["scc"])

                        if lrs[1] in lr_dists[lrs[0]]:
                            #curr_dist = np.mean(lr_dists[lrs[0]][lrs[1]][1])
                            #if abs(sr_dist - curr_dist) > 2000:
                            #    print("\t".join(["curr_dist: " + str(curr_dist), "sr_dist: " + str(sr_dist), "ctg1: " + ctg1, "ctg2: " + ctg2, lrs[0], lrs[1]]))
                            lr_dists[lrs[0]][lrs[1]][1].append(sr_dist)
                        else:
                            #print(lrs[0] + " + " + lrs[1] + ": " + ctg1 + " " + ctg2)
                            lr_dists[lrs[0]][lrs[1]]= ([],[sr_dist])
                        if lrs[0] in lr_dists[lrs[1]]:
                            lr_dists[lrs[1]][lrs[0]][1].append(-sr_dist)
                        else:
                            lr_dists[lrs[1]][lrs[0]]= ([],[-sr_dist])
    
    #save the distance matrix
    with open(dm_path, 'wb') as f:
        pickle.dump(lr_dists, f, pickle.HIGHEST_PROTOCOL)
        
else:
    with open(dm_path, 'rb') as f:
        lr_dists = pickle.load(f)
                        
    
for lrid,lrdists in lr_dists.items():
    for lrid2,dists in lrdists.items():
        if dists[0] == [] and len(dists[1]) > 1:
            pass
            #print("interesting: " + str(lrid) + " " + str(lrid2) + " " + str(dists[1]))
        if len(dists[0]) > 1  and len(dists[1]) > 1:
            #if abs(np.mean(dists[0]) - np.mean(dists[1])) > 50:
            #    print("\t".join([str(lrid),str(lrid2),str(np.mean(dists[0])),  str(np.mean(dists[1]))]))
            if abs(np.mean(dists[0]) - np.mean(dists[1])) > 500:
                print("\t".join([str(lrid),str(lrid2),str(np.mean(dists[0])),  str(np.mean(dists[1]))]))
                print("\t".join(["","",str(dists[0]),  str(dists[1])]))

pp = PdfPages("foo.pdf")
entries = 0
for lrid,lrdists in lr_dists.items():
    for lrid2,dists in lrdists.items():
        entries += 1
        if dists[1] != 0:
            nm = np.mean(dists[1])
            distances = [x-nm for x in dists[1]]
            plot1 = plt.figure()
            plt.hist(distances, range(-2000, 2000, 200))
            pp.savefig(plot1)
        if entries > 100:
            pp.close()
            sys.exit(0)


#for lrid,dists in lr_dists.items():
#    print("-"*200)
#    print(lrid)
#    print(dists)

        
devs = []
for dists in all_dists:
    stdev = (np.std(dists))
    if stdev < 50:
        pass
#        print(dists)
    devs.append(stdev)

# and now let's build a simple matrix
# taking an arbitrary but fixed ordering
lr_keys = lr_dists.keys()
matrix = []

for nr1, lr1 in enumerate(lr_keys):
    row = []
    for nr2, lr2 in enumerate(lr_keys):
        if lr1 in lr_dists:
            if lr2 in lr_dists[lr1]:
                row.append(1)
            else:
                row.append(0)
        else:
            row.append(0)

    matrix.append(row)

plt.imsave('connectedness.png', np.array(matrix).reshape(len(lr_keys),len(lr_keys)), cmap=cm.gray)
#plt.plot(np.array(matrix), cmap=cm.gray)
#plt.show()
print(lr_dists["0cf11d19-fcbf-4685-901f-32c4259eaf85"])

cluster = 0
unvisited_nodes = set(lr_keys)
while unvisited_nodes:
    cluster += 1
    print("Cluster " + str(cluster))
    start_node = unvisited_nodes.pop()
    current_cluster = [start_node]
    current_index = 0
    while(current_index != len(current_cluster)):
        for lr2 in lr_dists[current_cluster[current_index]]:
            if lr2 not in current_cluster and lr2 in unvisited_nodes:
                current_cluster.append(lr2)
                unvisited_nodes.remove(lr2)
        current_index += 1
    print(current_cluster)
    

# get pairs of overlapping reads
#plt.hist(devs,150)
#plt.yscale('log', nonposy='clip')
#plt.show()
