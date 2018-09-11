from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from itertools import combinations


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("--mincontigs", type=int, default=2,help="Minimum number of contigs on long read for the read to be considered")
parser.add_argument("--summaryfile", help="Contig Distance Summary file")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("linename", help="Name of cell line")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=25)

args = parser.parse_args()

# global data structures
contigs = {}
contig2scaf = {}
scaffolds = {}

# in clusters similar scaffolds will be stored, for manual curation
clusters = set()

if args.summaryfile:
    with open(args.summaryfile) as f:
        for line in f:
            sline = line.split()
            if sline[1] == "NA":
                continue
            if int(sline[2]) < args.mindepth:
                continue
            ctg1 = sline[0].split("_")[0].strip("+").strip("-")
            ctg2 = sline[0].split("_")[1].strip("+").strip("-")
            ori1 = sline[0].split("_")[0][0]
            ori2 = sline[0].split("_")[1][0]
            distance = float(sline[1])
            dist = int(distance)


if args.blacklistfile:
    blacklist = {}
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            blacklist[sline[0]] = sline[1]

for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)
    

print("Nr. of scaffolds: " + str(len(contigs)))

class Node:
    segments = {} # id1: (lc1, rc1)
    sizes_nc = [] # sizes of contigs that determine the average (that were Not Cut)
    
    def get_right_nb(self):
        pass
    
    def get_distance_table(self):
        pass
        
    
    def get_avg(self):
        return np.mean(self.sizes_nc)

    def __init__(self, contiginfo):
        #lc1 = 
        pass

class Edge:
    sizes = {} # no need for non cut edges
    def get_avg(self):
        avg = np.mean(self.sizes.values())
    def __init__():
        pass

# a scaffold is a connected graph of contig-nodes 
class Scaffold:
    distanceMatrix = []
    length = 0 # length is defined by the average distances and the average node sizes 
    nodes = set([])

    def update_distance_matrix(self):
        pass
    
    def compare(self, scaf2):
        #compare their distance matrices
        pass

    def to_graph_string(self):
        pass
        # give a graph string in order to visualize this ting
    
    def __init__(self):
        pass
    
    @classmethod
    def init_from_longread(cls,longread):
        newinst = cls()
        [lrid, lr] = longread
        for m in lr["maps"]:
            Node(m)
        return newinst
    

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
        if rid in lreads:
            lreads[rid]["maps"].append(data)
        else:
            lreads[rid] = {}
            lreads[rid]["length"] = int(lenr)
            lreads[rid]["maps"] = [data]


# get interesting reads
greadst = {}
for rid in lreads:
    counter = 0
    for item in lreads[rid]["maps"]:
        if item["contig"].endswith(args.linename) and item["contig"] != "1036QBL":
            counter +=1
            if counter >= args.mincontigs:
                greadst[rid] = lreads[rid]
                break

# turn reads around if necessary
for rid, lr in greadst.items():
    bw = 0
    fw = 0
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
                mapping["scr"] = lr["length"] - mapping["ecr"]
                mapping["ecr"] = lr["length"] - mapping["scr"]
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - mapping["scc"]
        
    # turn around and redefine wrong contigs
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.linename): 
            if mapping["strand"] == 1: #define a new contigname and turn it around
                mapping["contig"] = mapping["contig"] + "rc"
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - mapping["scc"]
                mapping["strand"] = 0
        

for rid in greadst:
    nscaff = Scaffold.init_from_longread((rid,lreads[rid]))
    scaffolds[id(nscaff)] = nscaff

scaffolds = {}


# cluster np-reads 
print("scaffolding long reads ....")

if args.summaryfile:
    print("adding short reads ....")
        


sys.exit(0)
