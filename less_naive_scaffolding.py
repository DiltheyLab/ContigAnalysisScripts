from argparse import ArgumentParser
from Bio import SeqIO
import pickle
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import random
from operator import itemgetter
from itertools import combinations, cycle, product
from collections import defaultdict, deque
import svgwrite
from scaffold import Scaffold, Longreads, LongReadSVG
import networkx as nx
from statistics import mean

parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
#parser.add_argument("--paf", help="Input is paf file", action="store_true", default = False)
parser.add_argument("--mincontigs", type=int, default=2,help="Minimum number of contigs on long read for the read to be considered")
parser.add_argument("--summaryfile", help="Contig Distance Summary file")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("--mergefile", help="File that contains merging information to retrieve the sequence of scaffolds after.")
parser.add_argument("--contigsmergefile", help="File with information on what contigs should be merged.")
parser.add_argument("--final_lrs", help="File to store final Longreads")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("SVG", help="Scaffolds are drawn to this SVG file")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=25)

args = parser.parse_args()

reads = {}
greads = {}
contigs = {}
allcontigs = {}
contig2scaffold = defaultdict(list)
cluster_counter = 1

srneighs = dict()

if args.mergefile:
    with open(args.mergefile, "w+") as mergef:
       pass 

ctgpos = {}
pos = 0
for read in SeqIO.parse(args.contigfile, "fasta"):
    pos += 1
    contigs[read.id] = len(read.seq)
    allcontigs[read.id] = len(read.seq)
    ctgpos[read.id] = pos


blacklist = defaultdict(list)
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist[sline[1]] = float(sline[2])
            else:
                blacklist[sline[0]].append(sline[1])
#print(blacklist)

def get_other_relpos(relpos):
    if relpos == "left":
        return "right"
    return "left"

def add_neighs(ctg1, ctg2, relpos, dist):
    if ctg1 in srneighs:
        if ctg2 not in [x[0] for x in srneighs[ctg1][relpos]]:
            srneighs[ctg1][relpos].append((ctg2,dist))
        if ctg2 in [x[0] for x in srneighs[ctg1][get_other_relpos(relpos)]]:
            pass
            #print("Conflict for " + ctg1 + " and " + ctg2)
            #print(srneighs[ctg1])
    else:
        srneighs[ctg1] = {"left": [], "right": []}
        srneighs[ctg1][relpos] = [(ctg2, dist)]

if args.summaryfile:
    with open(args.summaryfile) as f:
        for line in f:
            sline = line.split()
            #if sline[1] == "NA":
            #    continue
            if int(sline[2]) < args.mindepth:
                continue
            ctg1 = sline[0].split("_")[0].strip("+").strip("-")
            ctg2 = sline[0].split("_")[1].strip("+").strip("-")
            if ctg1 in blacklist or ctg2 in blacklist:
                continue
            ori1 = sline[0].split("_")[0][0]
            ori2 = sline[0].split("_")[1][0]
            distance = float(sline[3])
            dist = int(distance)
            if distance < 0 and abs(distance) > allcontigs[ctg1]:
                continue
            if distance < 1000 and distance > -1000:
                if ori1 == ori2:
                    if ori1 == "+":
                        add_neighs(ctg1,ctg2,"right",dist)
                        add_neighs(ctg2,ctg1,"left",dist)
                    else:
                        add_neighs(ctg2,ctg1,"right",dist)
                        add_neighs(ctg1,ctg2,"left",dist)
                else: # some contigs are defined in wrong orientation
                      # we still want to handle these
                    if ori1 =="+": # so ori2 == "-"  |ctg1 >| |< ctg2|
                        add_neighs(ctg1,ctg2,"right",dist)
                        add_neighs(ctg2,ctg1,"right",dist)
                    else: # |ctg1 <| |> ctg2|
                        add_neighs(ctg1,ctg2,"left",dist)
                        add_neighs(ctg2,ctg1,"left",dist)

# get rid of prefixes a_ b_ c_ and the like
def shortname(ctgname):
    if "_" in ctgname:
        return ctgname.split("_")[1]
    else:
        return ctgname



print("Nr. of Contigs: " + str(len(contigs)))



def cluster(scores, dists,scafs, verbose= False):
    # greedyly expand components with best matching thing
    gr = nx.DiGraph()
    if verbose:
        print(scores)
    for lr1i, lr1 in scores.items():
        if lr1i not in gr.nodes():
            gr.add_node(lr1i)
        ms = max(lr1.values())
        lr2i = max(lr1.keys(), key=lambda x: lr1[x])
        if verbose:
            print("lr1i: " + str(lr1i))
            print("lr2i: " + str(lr2i))
            print("ms: " + str(ms))
        if lr2i not in gr.nodes():
            gr.add_node(lr2i)
        #print(lr2i)
        #print(dists[lr1i][lr2i])
        if dists[lr1i][lr2i] is not None: # could be None if score was 0
            gr.add_edge(lr1i,lr2i,dist=dists[lr1i][lr2i])
            gr.add_edge(lr2i,lr1i,dist=dists[lr2i][lr1i])

    # order components
    creads = []
    for component in list(nx.connected_components((gr.to_undirected()))):
        # get all distances from random anchor node
        if len(component) == 1:
            creads.append(deque(component))
            continue
        anchor = random.sample(component,1)[0]
        worklist = deque([anchor])
        ns = gr.neighbors(anchor)
        for n in ns:
            worklist.append(n)
        while worklist:
            #print("worklist: " + str(len(worklist)))
            cn = worklist.popleft()
            ns = nx.neighbors(gr,cn)
            #print("ns: " + str(len(ns)))
            for n in ns:
                if (anchor, n) not in gr.edges():
                    gr.add_edge(anchor,n, dist=gr[anchor][cn]["dist"] + gr[cn][n]["dist"])
                    gr.add_edge(n, anchor, dist= - gr[anchor][cn]["dist"] - gr[cn][n]["dist"])
                    worklist.append(n)
        # now get all distances from the leftmost node
        minnode = min(gr[anchor], key = lambda x: gr[anchor][x]["dist"])
        maxnode = max(gr[anchor], key = lambda x: gr[anchor][x]["dist"]+scafs.lreads[x]["length"])
        print(minnode)
        print(maxnode)

        # let's jump from read to read 
        #reachable_nodes = gr[minnode]
        #while maxnode not in reachable_nodes:
        #mi
        # 
        
        offsets = scafs.get_possible_offsets(minnode,maxnode)
        pairscores = []
        for offset in offsets:
            pairscores.append(scafs.pseudoalign(minnode,maxnode,offset[1]-offset[0]))
        if offsets:
            print(max(pairscores))

        for node in set(component):
            gr.add_edge(minnode, node, dist = gr[minnode][anchor]["dist"] + gr[anchor][node]["dist"])
            gr.add_edge(node, minnode, dist = -gr[minnode][anchor]["dist"] - gr[anchor][node]["dist"])
        # put this in data structure that's used by the plotting, update dists
        creads.append(deque(sorted(component, key=lambda x: gr[minnode][x]["dist"])))
        #print(g[minnode])
    return [creads, gr]


# poor man's kmeans clustering
def split_distances(distances, tolerance):
    clusters = []
    indices = []
    for didx, dist in enumerate(distances):
        for cidx, cluster in enumerate(clusters):
            if abs(mean(cluster)-dist) < tolerance:
                cluster.append(dist)
                indices[cidx].append(didx)
                break
        else:
            clusters.append([dist])
            indices.append([didx])
    return (indices, clusters)

def merge_clusters(longreads, clustered_reads, distance_graph, minlength=1):
    # merge reads 
    pseudolongreads = dict()
    plr_nr = 1
    for cluster in clustered_reads:
        if len(cluster) < minlength:
            continue
        # perform Pseudo MSA and get consensus
        #print("Getting Pseudo-Longread for cluster...  (length: " + str(len(cluster)) + ")")
        pseudo_read_name = "pseudo_lr_" + str(plr_nr)
        plr_nr += 1
        pseudolongreads[pseudo_read_name] = dict()
        p = pseudolongreads[pseudo_read_name]
        p["lrids"] = cluster
        p["maps"] = []

        if len(cluster) == 1:
            p["maps"] = longreads.lreads[cluster[0]]["maps"].copy()
            #p["length"] = longreads.lreads[cluster[0]]["maps"].copy()
            continue
        
        origin = cluster[0]

        # pseudo MSA with consensus
        vstarts = defaultdict(list)
        ctgs = defaultdict(list)
        for lrid in cluster:
            lr = longreads.lreads[lrid]
            for contig in lr["maps"]:
                vstarts[contig["name"]].append(distance_graph[origin][lrid]["dist"] + contig["scr"] - contig["scc"])
                ctgs[contig["name"]].append(contig)

        for ctgn, ctgdists in vstarts.items():
            indices, clustered_dists = split_distances(ctgdists, 4000)
            #print(indices)
            for clusteridx, cluster in enumerate(clustered_dists):
                sccs = []
                eccs = []
                strands = []
                for ctgidx in indices[clusteridx]:
                    sccs.append(ctgs[ctgn][ctgidx]["scc"])
                    eccs.append(ctgs[ctgn][ctgidx]["ecc"])
                    strands.append(ctgs[ctgn][ctgidx]["strand"])
                newctg = {"strand": round(mean(strands)), "name": ctgn, "scc":round(mean(sccs)), "ecc":round(mean(eccs)), "scr":round(mean(cluster)) + round(mean(sccs)), "ecr":round(mean(cluster))+ round(mean(eccs))}
                p["maps"].append(newctg)
    return pseudolongreads

scafs = Longreads(args.inputfiles, blacklist, args.linename)
#scafs.filter_whitelist_ctgs(set(["1115APD"]))
scafs.filter_small_contigs(300)
scafs.filter_reverse_small_contigs(600)
scafs.filter_low_quality_contigs(0.81)
scafs.turn_longreads_around()
scafs.sort_by_starts()
scafs.filter_contigcounts(args.mincontigs)
#scafs.copy()

print("Nr. of reads: " + str(len(scafs.lreads)))

status = 0
for iteration in range(10):
    print("Pseudoaligning all... ", end="")
    lr_scores, lr_dists = scafs.pseudoalign_all()
    print("finished.")
    print("Clustering... ", end="")
    creads, gr = cluster(lr_scores, lr_dists, scafs)
    if len(scafs.lreads) == len(creads):
        status += 1
    print("finished.")
    print("Nr. of clusters: " + str(len(creads)))
    ps = merge_clusters(scafs, creads, gr, 1) if iteration != 0 else merge_clusters(scafs, creads, gr, 2)
    scafs = Longreads.init_from_dict(ps, scafs.cellline, contigs, scafs.lreads)
    scafs.sort_by_starts()
    
    if status == 2:
        scafs.filter_small_double_contigs(contigs, 0.85, False)
        scafs.filter_reverse_double_contigs(contigs, False)
        scafs.filter_overlapped_contigs(True)
    elif status == 4:
        scafs.filter_reverse_double_contigs(contigs, True)
        break

    for lrn in list(scafs.lreads.keys()):
        if not scafs.lreads[lrn]["length"]:
            scafs.delete(lrn)
                
        

lrsvg = LongReadSVG(args.SVG, zoom=200)
img = lrsvg.dwg
scafs.to_SVG(lrsvg, scafs.lreads.keys(), contigs, 30, lrsvg.yp, ctg_y_drawsize = 20)
img.save()

if args.final_lrs:
    with open(args.final_lrs, 'wb') as f:
        pickle.dump(scafs, f, pickle.HIGHEST_PROTOCOL)

sys.exit()
