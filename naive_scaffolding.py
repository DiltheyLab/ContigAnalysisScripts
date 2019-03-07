from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from itertools import combinations
from scaffold import Scaffolds
from collections import defaultdict


parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input files", nargs="+")
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("linename", help="Cell Line")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=20)
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")

args = parser.parse_args()

reads = {}
greads = {}
cgreads = []
contigs = {}


for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

print("Nr. of scaffolds: " + str(len(contigs)))

blacklist = defaultdict(list)
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist[sline[1]] = "y"
            else:
                blacklist[sline[0]].append(sline[1])

lrs = Scaffolds(args.inputfiles, blacklist, args.linename)
lrs.filter_contigcounts(2)
lrs.turn_longreads_around()
lrs.sort_contigs_in_reads()
greads = lrs.lreads



# cluster np-reads 
print("scaffolding long reads ....")
contig2cluster = {}
creads = {}
clusternr = 0
while len(greads) > 0:
    clusternr += 1
    current_cluster = {}
    current_contigs = set()
    # take a random read and build a cluster from it
    cr = greads.popitem()
    current_cluster[cr[0]] = cr[1]
    olen = 0
    while len(current_cluster) != olen:
        olen = len(current_cluster)
        for contig in cr[1]["maps"]:
            if not contig["name"].startswith("chr"):
                contigs.pop(contig["name"], None)
                current_contigs.add(contig["name"])
                contig2cluster[contig["name"]] = clusternr
        for readid,readval in greads.items():
            contig_found = False
            for contig in readval["maps"]:
                if not contig["name"].startswith("chr"):
                    if contig["name"] in current_contigs:
                        contig_found = True
                        current_cluster[readid] = readval
                        cr = (readid, greads.pop(readid))
                        break
                   
            if contig_found:
                break
    creads[clusternr] = current_cluster
print("Nr. of scaffolds: " + str(clusternr+len(contigs)) + " (" + str(clusternr) + " cluster + " + str(len(contigs))+ " contigs)")

# simpler data structure to collect contigs into scaffolds
scaffolds={}
for i, cluster in creads.items():
    current_contigs = set([])
    for readid,read in cluster.items():
        #print(read)
        for contig in read["maps"]:
            if not contig["name"].startswith("chr"):
                current_contigs.add(contig["name"])
    scaffolds[i] = current_contigs
#print(scaffolds)

def addcontig(ctg, cluster):
    contig2cluster[ctg] = cluster
    scaffolds[cluster].add(ctg)
    contigs.pop(ctg, None)
    
def mergecluster(cluster1, cluster2):
    for contig in scaffolds[cluster2]:
        contig2cluster[contig] = cluster1
        scaffolds[cluster1].add(contig)
    scaffolds.pop(cluster2)

# very lenient clustering of short reads
print("scaffolding short reads ....")
with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        if len(sline[0].split("_")) < 2:
            print("Problem with line: " + str(line.rstrip()))
        ctg1 = sline[0].split("_")[0].strip("+").strip("-")
        ctg2 = sline[0].split("_")[1].strip("+").strip("-")
        if sline[1] == "NA":
            continue
        if int(sline[2]) < args.mindepth:
            continue
        #moddist = float(sline[1])
        if ctg1 in contig2cluster:
            if ctg2 in contig2cluster:
                if contig2cluster[ctg1] != contig2cluster[ctg2]:
                    #print("merging clusters " + str(contig2cluster[ctg1]) + " and " + str(contig2cluster[ctg2]))
                    mergecluster(contig2cluster[ctg1], contig2cluster[ctg2])
            else:
                addcontig(ctg2, contig2cluster[ctg1])
            #print("cluster of ctg1 (" + ctg1 + "): " + str(contig2cluster[ctg1]))
        elif ctg2 in contig2cluster:
            addcontig(ctg1, contig2cluster[ctg2])
        else:
            clusternr += 1
            #print("new cluster: " + str(clusternr))
            scaffolds[clusternr] = set([ctg1, ctg2])
            contig2cluster[ctg1] = clusternr
            contig2cluster[ctg2] = clusternr
            contigs.pop(ctg1, None)
            contigs.pop(ctg2, None)
            
            
print("Nr. of scaffolds: " + str(len(scaffolds)+len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")
#print(contigs)
for scaf in scaffolds:
    print(scaffolds[scaf])
    print("-"*200)

def get_rightmost_contig(scaf_id):
    print(creads[scaf_id])
    if len(creads[scaf_id]) > 1:
        print("not implemented yet")
    else:
        rid, read = creads[scaf_id].items()
        for ov in read["maps"]:
            print(ov)

#get_rightmost_contig(1)
#sys.exit(0)
