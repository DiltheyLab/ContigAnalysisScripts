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

class Scaffold:
    
    @classmethod
    def init_from_longread(cls,lr):
        newinst = cls()
        return newinst
    

# nanopore reads
lreads = {}
with open(args.efile) as f:
    for line in f:
        #sline = line.split()
        [rid, ctg, t2, t3, t4, scr, ecr, lenr, strand, scc, ecc, lenc, t12, t13, t14, t15, t16] = line.split()
        payload = {"contig":ctg,"strand":strand,"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
        if rid in blacklist:
            if blacklist[rid] == ctg:
                continue
        if rid in lreads:
            lreads[rid]["maps"].append(payload)
        else:
            lreads[rid] = {}
            lreads[rid]["length"] = int(lenr)
            lreads[rid]["maps"] = [payload]

scaffolds = {}
# get interesting reads
# and sort contigs by left coordinate
greadst = {}
for rid in lreads:
    counter = 0
    for item in lreads[rid]["maps"]:
        if item["contig"].endswith(args.linename) and item["contig"] != "1036QBL":
            counter +=1
            if counter >= args.mincontigs:
                greadst[rid] = reads[rid]
                break

for rid in greadst:
    nscaff = Scaffold.init_from_longread((rid,reads[rid]))
    scaffolds[id(nscaff)] = nscaff


# cluster np-reads 
print("scaffolding long reads ....")

if args.summaryfile:
    print("adding short reads ....")
        


sys.exit(0)
