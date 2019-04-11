from argparse import ArgumentParser
from Bio import SeqIO
import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import combinations, cycle, product
from collections import defaultdict, deque
from scaffold import Scaffold, Longreads, LongReadSVG
from statistics import mean

parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("sequencefile", help="Input File in FASTQ format")
parser.add_argument("contigstringfile", help"Input File containing subsequent longread ids")
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

firstcontig = ""
with open(args.contigstringfile) as f:
    for line in f:
        if len(line.split()) != 2:
            continue
        else:
            firstcontig = line.split()[1].rstrip() + args.linename
            break
    else:
        print("Error while reading " + args.contigstringfile)
        print("Could not find first contig")
        sys.exit()
            


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
            items.append([lrid, contig])
        elif len(len.split()) == 1:
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


out_sequence = contigs[firstcontig]

for item1,item2 in zip(items[:-1],items[1:]):
    lastcontig = contig + args.linename
    status = 0
    for ctg in scafs.lreads[lrid]["maps"]:
        if status == 0:
            if ctg["name"] == firstcontig:
                status = 1
                len(contigs[c

            else:
                continue
        elif status == 1:
            if ctg["name"] == lastcontig
        
        
    elif len(line.split()) == 1:
        ctg = line + args.linename
        
    
sys.exit()
    

