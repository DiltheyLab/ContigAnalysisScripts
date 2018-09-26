from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from itertools import combinations
import svgwrite


parser = ArgumentParser()
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("--mindepth", help="Minimal depth in the short read summary file.", type=int, default=25)

args = parser.parse_args()

srneighs = dict()

def get_other_relpos(relpos):
    if relpos == "left":
        return "right"
    return "left"
def add_neighs(ctg1, ctg2, relpos, dist):
    if ctg1 in srneighs:
        if ctg2 not in [x[0] for x in srneighs[ctg1][relpos]]:
            srneighs[ctg1][relpos].append((ctg2,dist))
        if ctg2 in [x[0] for x in srneighs[ctg1][get_other_relpos(relpos)]]:
            print("Conflict for " + ctg1 + " and " + ctg2)
            print(srneighs[ctg1])
    else:
        srneighs[ctg1] = {"left": [], "right": []}
        srneighs[ctg1][relpos] = [(ctg2, dist)]

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
        if distance < 1000 and distance > -200:
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

contigs = {}

for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

print("Nr. of scaffolds: " + str(len(contigs)))

reads = {}

# nanopore read info
with open(args.efile) as f:
    for line in f:
        sline = line.split()
        rid = sline[0]
        ctg = sline[1]
        strand = int(sline[8])
        scr = int(sline[5])
        ecr = int(sline[6])
        lenr = int(sline[7])
        scc = int(sline[9])
        ecc = int(sline[10])
        lenc = int(sline[11])
        payload = {"contig":ctg,"strand":strand,"scr":scr,"ecr":ecr,"scc":scc,"ecc":ecc,"lenc":lenc}
        if rid in reads:
            reads[rid]["maps"].append(payload)
        else:
            reads[rid] = {}
            reads[rid]["length"] = int(sline[7])
            reads[rid]["maps"] = [payload]

scaffolds = {}
greadst = {}
greads = {}

for rid in reads:
    counter = 0
    for item in reads[rid]["maps"]:
        if item["contig"].endswith(args.linename):
            counter +=1
            if counter >= 2:
                greadst[rid] = reads[rid]
                break

# sort contigs by left coordinate
for rid in greadst:
    #print(reads[rid]["overlaps"])
    soverlaps = sorted(reads[rid]["maps"], key = itemgetter("scr"))
    greads[rid] = greadst[rid]
    greads[rid]["maps"]=soverlaps

lrneighs = {}

def add_lr_neighs(ctg1, ctg2, relpos, dist):
    if ctg1 in lrneighs:
        if lrneighs[ctg1][relpos] != []:
            current_neighbour = lrneighs[ctg1][relpos][0] 
            if current_neighbour[1] > dist:
                lrneighs[ctg1][relpos] = [(ctg2, dist)]
        else:
            lrneighs[ctg1][relpos]=[(ctg2,dist)]
    else:
        lrneighs[ctg1] = {"left": [], "right": []}
        lrneighs[ctg1][relpos] = [(ctg2, dist)]


for rid, read in greads.items():
    for part1,part2 in zip(read["maps"], read["maps"][1:]):
        ctg1, ctg2 = [part1["contig"],part2["contig"]]
        if ctg1.endswith(args.linename) and ctg2.endswith(args.linename):
            distance = part2["scr"] - part1["ecr"]- (part1["lenc"] - part1["ecc"]) - part2["scc"] + 1
            if part1["strand"] == 0:
                add_lr_neighs(ctg1,ctg2,"right",distance)
            else:
                add_lr_neighs(ctg1,ctg2,"left",distance)
            if part2["strand"] == 0:
                add_lr_neighs(ctg2,ctg1,"left",distance)
            else:
                add_lr_neighs(ctg2,ctg1,"right",distance)



# we have srneighs and lrneighs
#print(srneighs)
l = "left"
r = "right"

for ctg,neighs in lrneighs.items():
    if len(neighs[l])> 1:
        print("contig " + ctg + " has left neighbours " + str(neighs[l]))
        continue
    if len(neighs[r])> 1:
        print("contig " + ctg + " has right neighbours " + str(neighs[r]))
        continue
    
    if neighs[l] != []:
        leftn, leftd = neighs[l][0]
    else:
        continue
#    if neighs[r] != []:
#        rightn = neighs[r][0][0]
    if ctg in srneighs and leftn not in  [x[0] for x in srneighs[ctg][l]] :
        print(str(ctg) + " has no " + str(leftn) + " in the short reads. " + str(leftd) + " distance")

