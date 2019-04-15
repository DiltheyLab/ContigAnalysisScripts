from argparse import ArgumentParser
from Bio import SeqIO
import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import combinations, cycle, product
from collections import defaultdict, deque
from scaffold import Scaffold, Longreads, LongReadSVG, revcomp
from statistics import mean

parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("sequencefile", help="Input File in FASTQ format")
parser.add_argument("contigstringfile", help="Input File containing subsequent longread ids")
parser.add_argument("contigfile", help="Contig File in FASTA format")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("outfile", help="Sequence Output File")
args = parser.parse_args()

contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = str(read.seq)

lreadseqs = {}
for read in SeqIO.parse(args.sequencefile, "fastq"):
    lreadseqs[read.id] = str(read.seq)

print("Nr. of Contigs: " + str(len(contigs)))

scafs = Longreads(args.inputfiles, None, args.linename)
#scafs.filter_small_contigs(300)
#scafs.filter_reverse_small_contigs(600)
scafs.filter_low_quality_contigs(0.81)
scafs.turn_longreads_around()
scafs.sort_by_starts()

print("Nr. of reads: " + str(len(scafs.lreads)))

#with open(args.contigstringfile) as f:
#    for line in f:
#        if len(line.split()) != 2:
#            continue
#        else:
#            firstcontig = line.split()[1].rstrip() + args.linename
#            break
#    else:
#        print("Error while reading " + args.contigstringfile)
#        print("Could not find first contig")
#        sys.exit()
            


items = []
with open(args.contigstringfile) as f:
    for line in f:
        if len(line.split()) == 2:
            lrid, contig = line.split()
            if lrid not in scafs.lreads:
                print("Error while reading " + args.contigstringfile)
                print("Did not load " + str(lrid) + " from the input files")
                sys.exit()
            if lrid not in lreadseqs:
                print("Error while reading " + args.contigstringfile)
                print("Did not load sequence for " + str(lrid))
                sys.exit()
            ctgn = contig + args.linename
            if ctgn not in contigs:
                print("Error while reading " + args.contigstringfile)
                print("Contig " + str(contig) + " not found.")
                sys.exit()
            items.append([lrid, ctgn])
        elif len(line.split()) == 1:
            ctgn = line.rstrip() + args.linename
            if ctgn not in contigs:
                print("Error while reading " + args.contigstringfile)
                print("Contig " + str(contig) + " not found.")
                sys.exit()
            items.append([ctgn])

        else:
            print("Error while reading " + args.contigstringfile)
            print("The following line does not have the appopriate format")
            print(line)


for item in items:
    if len(item) == 2:
        if item[0] not in lreadseqs:
            print(str(item[0]) +" not found in sequence file. Aborting !")
            sys.exit()
distances = {}
distances[("2000APD","1503APD")] = -20
distances[("1503APD","476APD")] = -13
distances[("476APD","1405APD")] = -9
distances[("2316APD","1498APD")] = -12
distances[("504APD","49APD")] = -447
distances[("849APD","1471APD")] = -11
firstcontig = "2345APD"

out_sequence = contigs[firstcontig]
last_ctgn = firstcontig
print("Getting sequence ...")

for item in items:
    #print("length of sequence: " + str(len(out_sequence)))
    #print(item)
    if len(item) == 2:
        first_ctgn = last_ctgn
        lrid, last_ctgn = item
        if "reverse" in scafs.lreads[lrid]:
            #print(lrid + " is reverse complimentary")
            lr_seq = revcomp(lreadseqs[lrid])
        else:
            lr_seq = lreadseqs[lrid]
        status = 0
        for ctg in scafs.lreads[lrid]["maps"]: # they should be ordered
            #print("ctgn: " + ctg["name"])
            if ctg["strand"] == 1: # no reverse contigs allowed
                continue
            if status == 0:
                if ctg["name"] == first_ctgn:
                    status = 1
                    last_used_ctg = ctg.copy()
                    # get things
                else:
                    continue
            elif status == 1:
                tocut = len(contigs[last_used_ctg["name"]]) - last_used_ctg["ecc"]
                #print("to cut: " +str(tocut) + " " + last_used_ctg["name"])
                if tocut > 0:
                    out_sequence = out_sequence[:-tocut]
                if last_used_ctg["ecr"] > ctg["scr"]:
                    out_sequence = out_sequence[:ctg["scr"]-last_used_ctg["ecr"]]
                    out_sequence += contigs[ctg["name"]][ctg["scc"]-1:]
                else:
                    out_sequence += lr_seq[last_used_ctg["ecr"]-1:ctg["scr"]-1]
                    out_sequence += contigs[ctg["name"]][ctg["scc"]-1:]
                if ctg["name"] == last_ctgn:
                    break
                last_used_ctg = ctg
    elif len(item) == 1:
        if last_ctgn == item[0]:
            continue
        else:
        # get distance between last_ctg and item[0]
            if (last_ctgn,item[0]) not in distances:
                print("get distance between " + last_ctgn + " and " + item[0])
                sys.exit()
            else:
                out_sequence = out_sequence[:distances[(last_ctgn, item[0])]]
                out_sequence += contigs[item[0]]
                last_ctgn = item[0]
                
            
with open(args.outfile, "w") as out:
    out.write(">MHC_APD\n")
    out.write(out_sequence)

        
    

