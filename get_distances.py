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
parser.add_argument("summaryfile", help="Contig distance summary file")
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
            
            pairs[(get_full_name(ctg1),get_full_name(ctg2))] = {"lr":[], "id":[], "sr":[]}

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
    return sorted(distances, key = lambda x: abs(x))[0] # returns the value that is |abs|-smallest

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
    print(str(pair) + "\t" + str(pairs[pair]["lr"]))

sys.exit()

# get distances from short reads
with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        ctg1, ctg2 = sline[0].split("_")
        if ctg1.startswith("-") or ctg2.startswith("-"):
            continue
        if sline[1] == "NA":
            continue
        #if float(sline[4]) > 2:
        #    continue
        moddist = float(sline[1])
        #if int(ctg1.rstrip(args.cellline)) < int(ctg2.rstrip(args.cellline)):
        pair = (ctg1.strip("+").strip("-"), ctg2.strip("+").strip("-"))
        if pair in pairs:
            pairs[pair]["sr"].append(moddist)
        #else:
        #    pairs[pair] = [moddist]
#pairs[("224APD","976APD")] = {"lr": [], "sr": [16.67]}
#pairs[("916APD","1247APD")] = {"lr": [], "sr": [116.8]}
#pairs[("2406APD","1671APD")] ={"lr": [], "sr": [52.3]}
#pairs[("105APD","726APD")] ={"lr": [], "sr": [-40.1]}
#pairs[("726APD","1674APD")] ={"lr": [], "sr": [304.0]}
#pairs[("2080APD","2377APD")] ={"lr": [], "sr": [361.9]}
#pairs[("2377APD","928APD")] ={"lr": [], "sr": [75.6]}

# weird things from another planet
print("Problem pairs")
for pair in pairs:
    if pairs[pair]["sr"] == [] and pairs[pair]["lr"] == []:
        print(pair)
print("-"*50)

# get alignment scores
scoret = pickle.load(open(args.alignment_scores, "rb"))
scores = {}
for pair in scoret: 
    scores[pair] = [float(x) for x in scoret[pair]]
    

lzpairs = {}
for pairn, pair in pairs.items():
    if len(pair["lr"]) > 0:
        lrm = median(pair["lr"])
    else:
        lrm = None
    if len(pair["sr"]) > 0:
        srm = sum(pair["sr"])/len(pair["sr"])
    else:
        srm = None
    #if (not lrm  or lrm <= 0) and (not srm or srm <= 0) and (lrm or srm):
    pairs[pairn]["lrm"] = lrm
    pairs[pairn]["srm"] = srm
    if (srm and srm <= 0) or (lrm and lrm<= 0):
        lzpairs[pairn] = srm

def get_next_good_distance(pairn, dist, margin):
    dlen = len(scores[pairn])
    if pairn not in scores:
        return None
    maxv = 0
    maxi = None
    for i in range(margin):
        if dist+i+dlen < dlen and dist+i+dlen >= 0:
            if scores[pairn][dist+i+dlen] > 0.95 and scores[pairn][dist+i+dlen] > maxv:
                maxv = scores[pairn][dist+i+dlen]
                maxi = dist + i
        if dist-i+dlen < dlen and dist-i+dlen >= 0:
            if scores[pairn][dist-i+dlen] > 0.95 and scores[pairn][dist-i+dlen] > maxv:
                maxv = scores[pairn][dist-i+dlen]
                maxi = dist - i
    return maxi

def search_forward(pairn, dist):
    dlen = len(scores[pairn])
    if pairn not in scores:
        return None
    maxi = 0
    maxv = 0
    for i in range(dist,0):
        if scores[pairn][i+dlen] > 0.95 and scores[pairn][i+dlen] > maxv:
            maxv = scores[pairn][dist+i+dlen]
            maxi = dist + i
    return maxi


real_distances = {}
for pairn, pair in pairs.items():
    if pairn not in lzpairs:
        continue
    if pair["lrm"] and pair["lrm"] <= 0:
        d = get_next_good_distance(pairn, int(round(pair["lrm"])), 40) 
        if d:
            real_distances[pairn] = d
    if pair["srm"] and pair["srm"] <= 0 and pair["lr"] == []:
        d = get_next_good_distance(pairn, int(round(pair["srm"])), 40) 
        if d:
            real_distances[pairn] = d
    if pairn not in real_distances:
        for dpred in [int(x) for x in pair["lr"]]:
            if dpred < 0:
                d = get_next_good_distance(pairn, dpred, 40) 
                if d:
                    real_distances[pairn] = d
                    break
        if pair["lr"] == [] and pair["srm"] < 0:
            d = get_next_good_distance(pairn, int(round(pair["srm"])), 40) 
            if d:
                real_distances[pairn] = d
            else:
                real_distances[pairn] = int(round(pair["srm"]))
    if pairn not in real_distances and pair["lrm"] != None and mean(pair["lr"]) > -200:
        real_distances[pairn] = search_forward(pairn, int(round(pair["lrm"])))
    elif pairn not in real_distances and pair["lrm"] != None:
        #TODO use the contig mapping information from LRs to take the "correct" sequence
        real_distances[pairn] = int(round(pair["lrm"]))
        

print("total pairs: " + str(len(pairs)))
print("lower zero pairs: " + str(len(lzpairs)))
print("lower zero pairs solved: " + str(len(real_distances)))
for pair in lzpairs.keys() - real_distances.keys():
    print(str(pair) + ": " + str(pairs[pair]))
#print(real_distances)
#print("-"*30)
#with open(args.distancefile,"w+") as out:
#    for pairn,paird in real_distances.items():
#        out.write("\t".join([pairn[0], pairn[1], str(paird)]) + "\n")
#    for pair in set(lzpairs.keys() - real_distances.keys()):
#        out.write("\t".join([pair[0], pair[1], str(pairs[pair])]) + "\n")
#    for pair in set(pairs.keys() - real_distances.keys()):
#        out.write("\t".join([pair[0], pair[1], str(pairs[pair])]) + "\n")


for pair in set(pairs.keys() - real_distances.keys()): # all pairs with positive distance
    ctg1, ctg2 = pair
    ids1 = contig2lrid[ctg1]
    ids2 = contig2lrid[ctg2]
    if pairs[pair]["lr"] == []:
        #print(pairs[pair]["sr"])
        real_distances[pair] = pairs[pair]["srm"]

print("lz and non-resolvables: " + str(len(real_distances)))

###############################################################
# Okay now let's get the sequences for the positive distances #
###############################################################
pathtocontigs = "/home/houwaart/Projects/ImmunoPore/APD/APDContigs.fasta"
pathtofastq = "/home/houwaart/Projects/ImmunoPore/APD/155-NGS-090818_both.mapsToContigs.fastq"

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
            if lrid in reads:
                if reads[lrid]["turned"]:
                    reads[lrid]["seq"] = revcomp(line.rstrip())
                else:
                    reads[lrid]["seq"] = line.rstrip()
print("done.")


for pair in set(pairs.keys() - real_distances.keys()): # all pairs with positive distance
    ctg1, ctg2 = pair
    #ids1 = contig2lrid[ctg1]
    #ids2 = contig2lrid[ctg2]
    print(pair , ": " ,pairs[pair], end="\r")
    #print(ids1&ids2)
    if pairs[pair]["lr"] == []:
        #print(pairs[pair]["sr"])
        continue
    with open("infixes.fasta", "w") as f:
        f.write("")
    for rid,dist in zip(pairs[pair]["id"], pairs[pair]["lr"]):
        if dist < 0:
            continue
        for contig in reads[rid]["overlaps"]:
            if ctg1 == contig["contig"]:
                pleft = contig["ecr"] - 100
                left = pleft if pleft > 0 else 0
            if ctg2 == contig["contig"]:
                pright = contig["scr"] + 100
                right = pright if pright < reads[rid]["length"] else reads[rid]["length"]
        with open("infixes.fasta", "a") as f:
            f.write(">" + rid +"_infix\n")
            f.write(reads[rid]["seq"][left:right] + "\n")
    
    # write 
    with open("tmp1.fasta", "w") as f:
        f.write(">" + ctg1 +"\n")
        f.write(contigs[shortname(ctg1)] + "\n")
    with open("tmp2.fasta", "w") as f:
        f.write(">" + ctg2 +"\n")
        f.write(contigs[shortname(ctg2)] + "\n")
    points = subprocess.run(["./alignments3","infixes.fasta", "tmp1.fasta", "tmp2.fasta"], stdout=subprocess.PIPE, universal_newlines=True)
    #print(points.stdout)
    results = points.stdout.split("\n")
    consensus = results[-2]
    real_distances[pair] = consensus
    #sys.exit(0)
print(len(real_distances))


with open(args.distancefile,"wb") as out:
    pickle.dump(real_distances, out, pickle.HIGHEST_PROTOCOL)

