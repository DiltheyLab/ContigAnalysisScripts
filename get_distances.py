from argparse import ArgumentParser
from operator import itemgetter
from collections import defaultdict, OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import sys
import pickle
from itertools import combinations, product
from statistics import median, mean
from Bio import SeqIO
from scaffold import Scaffolds


parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("pairfile", help="File of contig pairs")
parser.add_argument("real_distances", help="Contig distances file")
parser.add_argument("alignment_scores", help="File containing alignment scores")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("distancefile", help="Output file of pairs with distances")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
args = parser.parse_args()

reads = {}
cgreads = []

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
lrs.filter_small_contigs(300)
#lrs.filter_reverse_small_contigs(600)

reverse_mappers = set()
reverse_mappers.add("344DBB")
reverse_mappers.add("472DBB")
lrs.turn_longreads_around(reverse_mappers)
lrs.sort_by_starts()

contig2lrid = lrs.ctg2lreads

def get_full_name(short_ctgn):
    if "_" in short_ctgn:
        nr = short_ctgn.split("_")[0]
        nr += "DBB_rc"
    else:
        nr = short_ctgn + "DBB"
    return nr

def sname(ctgn):
    if "_" in ctgn:
        return ctgn.split("_")[0]
    else:
        return ctgn


# get all pairs from file
pairs = OrderedDict()
with open(args.pairfile) as f:
    for line in f:
        if not line.startswith(">"):
            ctg1, ctg2 = line.rstrip().split()
            
            pairs[(get_full_name(ctg1),get_full_name(ctg2))] = {"lr":[], "id":[], "reald":[]}

def is_rc(ctgn):
    if ctgn.endswith("_rc"):
        return True
    else:
        return False

def pair_exists(read, ctgn1, ctgn2):
    for ctg1 in read["mapsc"][sname(ctgn1)]:
        if (is_rc(ctgn1) and ctg1["strand"] == 1) or (not is_rc(ctgn1) and ctg1["strand"] == 0):
            break
    else:
        return False
    for ctg2 in reversed(read["mapsc"][sname(ctgn2)]):
        if (is_rc(ctgn2) and ctg2["strand"] == 1) or (not is_rc(ctgn2) and ctg2["strand"] == 0):
            break
    else:
        return False

    if (ctg1["scr"] - ctg1["scc"]) < (ctg2["scr"] - ctg2["scc"]):
        return True
    else:
        return False
       
def get_smallest_distance(read, ctgn1, ctgn2):
    distances = []
    ec1rs = []
    sc2rs = []
    for ctg1, ctg2 in product(read["mapsc"][sname(ctgn1)],read["mapsc"][sname(ctgn2)]):
        if (is_rc(ctgn1) and ctg1["strand"] == 0) or (not is_rc(ctgn1) and ctg1["strand"] == 1):
            continue
        if (is_rc(ctgn2) and ctg2["strand"] == 0) or (not is_rc(ctgn2) and ctg2["strand"] == 1):
            continue
        start_rest1 = ctg1["scc"]
        start_rest2 = ctg2["scc"]
        end_rest1 = lrs.contig_lengths[ctg1["name"]] - ctg1["ecc"]
        end_rest2 = lrs.contig_lengths[ctg2["name"]] - ctg2["ecc"]
        vend1 = ctg1["ecr"] + end_rest1 if ctg1["strand"] == 0 else ctg1["ecr"] + start_rest1
        vbeg2 = ctg2["scr"] - start_rest2 if ctg2["strand"] == 0 else ctg2["scr"] - end_rest1
        distances.append(vbeg2 - vend1)
        ec1rs.append(ctg1["ecr"])
        sc2rs.append(ctg2["scr"])
    absdists = [abs(x) for x in distances]
    sidx = absdists.index(min(absdists))
    return (distances[sidx], ec1rs[sidx], sc2rs[sidx])
    #dist = sorted(distances, key = lambda x: abs(x))[0] # returns the value that is |abs|-smallest


def get_all_lrs_with_sig_from_pairs(pairs, lreads):
    good_reads = {}
    for rid, read in lreads.items():
        for ctg1, ctg2 in pairs:
            if sname(ctg1) in read["mapsc"] and sname(ctg2) in read["mapsc"]:
                if pair_exists(read, ctg1, ctg2):
                    good_reads[rid] = read
    return good_reads
            
            
good_reads = get_all_lrs_with_sig_from_pairs(pairs, lrs.lreads)

# get distances from long reads
for rid, read in good_reads.items():
    for pair in pairs:
        ctg1, ctg2 = pair
        if sname(ctg1) in read["mapsc"] and sname(ctg2) in read["mapsc"] and pair_exists(read,ctg1, ctg2):
            pairs[pair]["lr"].append(get_smallest_distance(read, ctg1, ctg2))
            pairs[pair]["id"].append(rid)

for pair in pairs:
    print(str(pair) + "\t" + str(pairs[pair]["lr"]) + "\t" + str(median(pairs[pair]["lr"][0])))


# long read we need to take care of
to_align = set()
real_distances = {}
with open(args.real_distances) as f:
    for line in f:
        ctg1, ctg2, dist = line.split()
        ctg1n = get_full_name(ctg1)
        ctg2n = get_full_name(ctg2)
        if ":" in dist:
            to_align.add((ctg1n,ctg2n,dist))
        else:
            pairs[(ctg1n,ctg2n)]["reald"] = int(dist)
            real_distances[(ctg1n,ctg2n)] = int(dist)

#print(pairs)

###############################################################
# Okay now let's get the sequences for the positive distances #
###############################################################
pathtocontigs = "/home/houwaart/Projects/ImmunoPore/DBB/DBBContigs.fasta"
pathtofastq = "/home/houwaart/Projects/ImmunoPore/DBB/dbb_lib1_lib2_lib3.fastq"

def shortname(ctg):
    if "_" in ctg:
        return ctg.split("_")[1]
    else:
        return ctg

print("loading contigs")
contigs = {}
for read in SeqIO.parse(pathtocontigs, "fasta"):
    contigs[read.id] = str(read.seq)
print("done.")


complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
def revcomp(instring):
    outstring = ""
    for char in instring:
        outstring += complement[char]
    return outstring[::-1]

print("Loading sequences ... ") 
with open(pathtofastq) as f:
    for i, line in enumerate(f):
        if i % 4 == 0:
            lrid = line.split(" ")[0][1:]
        elif i % 4 == 1:
            if lrid in lrs.lreads:
                if lrs.lreads[lrid]["turned"]:
                    lrs.lreads[lrid]["seq"] = revcomp(line.rstrip())
                else:
                    lrs.lreads[lrid]["seq"] = line.rstrip()
print("done.")


for item in to_align:
    ctg1, ctg2, dists = item
    #ids1 = contig2lrid[ctg1]
    #ids2 = contig2lrid[ctg2]
    pair = (ctg1,ctg2)
    #print(pair , ": " ,pairs[pair]["reald"], end="\r")
    dist_low , dist_high = [int(x) for x in dists.split(":")]
    with open("infixes.fasta", "w") as f:
        pass

    #print("-"*20)
    for distitem, rid in zip(pairs[pair]["lr"],pairs[pair]["id"]):
        dist, ec1r, sc2r = distitem
        if dist >= dist_low and dist <= dist_high:
            pleft = ec1r -100
            left = pleft if pleft > 0 else 0
            pright = sc2r + 100
            right = pright if pright < lrs.lreads[rid]["length"] else lrs.lreads[rid]["length"]
            print(dist, rid,left,right)

            with open("infixes.fasta", "a") as f:
                f.write(">" + rid +"_infix\n")
                f.write(lrs.lreads[rid]["seq"][left:right] + "\n")

    
    # write 
    with open("tmp1.fasta", "w") as f:
        f.write(">" + ctg1 +"\n")
        if ctg1.endswith("_rc"):
            f.write(revcomp(contigs[sname(ctg1)]) + "\n")
        else:
            f.write(contigs[sname(ctg1)] + "\n")
    with open("tmp2.fasta", "w") as f:
        f.write(">" + ctg2 +"\n")
        if ctg2.endswith("_rc"):
            f.write(revcomp(contigs[sname(ctg2)]) + "\n")
        else:
            f.write(contigs[sname(ctg2)] + "\n")
    points = subprocess.run(["./alignments3","infixes.fasta", "tmp1.fasta", "tmp2.fasta"], stdout=subprocess.PIPE, universal_newlines=True)
    #print(points.stdout)
    results = points.stdout.split("\n")
    consensus = results[-2]
    real_distances[pair] = consensus
    print(pair)
    print(consensus)
    #sys.exit(0)
print(len(real_distances))

with open(args.distancefile,"wb") as out:
    pickle.dump(real_distances, out, pickle.HIGHEST_PROTOCOL)
