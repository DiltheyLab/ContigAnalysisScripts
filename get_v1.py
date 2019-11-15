from argparse import ArgumentParser
from Bio import SeqIO
import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np
import random
from itertools import combinations, cycle, product
from collections import defaultdict, deque
from scaffold import Longreads, LongReadSVG, revcomp
from statistics import mean

parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("sequencefile", help="Input File in FASTQ format")
parser.add_argument("contigstringfile", help="Input File containing subsequent longread ids")
parser.add_argument("contigfile", help="Contig File in FASTA format")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("--dryrun", default=False, action="store_true", help="Do not output sequence to file, but only print the length of the resulting sequence to stdout.")
parser.add_argument("--logfile", help="Logging File")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("outfile", help="Sequence Output File")
args = parser.parse_args()

seq1="CCCTTTTTGGGAAA"
print(seq1)
print(revcomp(seq1))
contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = str(read.seq)

blacklist = defaultdict(list)
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist[sline[1]] = float(sline[2])
            else:
                blacklist[sline[0]].append(sline[1])

if not args.dryrun:
    lreadseqs = {}
    for read in SeqIO.parse(args.sequencefile, "fastq"):
        lreadseqs[read.id] = str(read.seq)

loghandle = None
if args.logfile:
    loghandle = open(args.logfile,"w+")

print("Nr. of Contigs: " + str(len(contigs)))

scafs = Longreads(args.inputfiles, blacklist, args.linename)
#scafs.filter_small_contigs(300)
#scafs.filter_reverse_small_contigs(600)
scafs.filter_low_quality_contigs(0.76)
#print(scafs.lreads["88e9f524-98be-4630-98bd-06054fd07355"]["maps"][10])
scafs.turn_longreads_around(logging=loghandle)
#print(scafs.lreads["88e9f524-98be-4630-98bd-06054fd07355"]["maps"][10])
#print(scafs.lreads["88e9f524-98be-4630-98bd-06054fd07355"])
scafs.sort_by_starts()
scafs.filter_contigs_by_coverage(0.5,ignore_ends=True, verbose = False)

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
            
seqlength = 0


items = []
with open(args.contigstringfile) as f:
    for line in f:
        if line.startswith("#"):
            continue
        if len(line.split()) == 2:
        
            lrid, contig = line.split()
            if lrid not in scafs.lreads:
                print("Error while reading " + args.contigstringfile)
                print("Did not load " + str(lrid) + " from the input files")
                sys.exit()
            if not args.dryrun:
                if lrid not in lreadseqs:
                    print("Error while reading " + args.contigstringfile)
                    print("Did not load sequence for " + str(lrid))
                    sys.exit()
            #ctgn = contig + args.linename
            ctgn = contig
            if ctgn not in contigs:
                print("Error while reading " + args.contigstringfile)
                print("Contig " + str(contig) + " not found.")
                sys.exit()
            items.append([lrid, ctgn])
        elif len(line.split()) == 1:
            #ctgn = line.rstrip() + args.linename
            ctgn = line.rstrip() 
            if ctgn not in contigs:
                print("Error while reading " + args.contigstringfile)
                print("Contig " + str(ctgn) + " not found.")
                sys.exit()
            items.append([ctgn])

        else:
            print("Error while reading " + args.contigstringfile)
            print("The following line does not have the appopriate format")
            print(line)


for item in items:
    if args.dryrun:
        break
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
distances[("1192DBB","462DBB")] = -9712
#firstcontig = "2345APD"
#firstcontig = "901DBB"
#firstcontig = "1003SSTO"
firstcontig = items[0][0]
print("First contig: " + firstcontig)
items = items[1:]

out_sequence = contigs[firstcontig]
overall_length = len(out_sequence)
last_ctgn = firstcontig

# check that important contigs weren't filtered out
for item1, item2 in zip(items[:-1], items[1:]):
    if len(item2) == 1:
        continue
    ctg1n = item1[0] if len(item1) == 1 else item1[1]
    lrid2, ctg2n = item2
    for ctgn in [ctg1n,ctg2n]:
        if ctgn not in scafs.lreads[lrid2]["ctgset"]:
            print("Problem found.")
            print(ctgn + " not in " + lrid2)
            sys.exit()
if args.dryrun:
    print("Getting length of sequence ...")
else:
    print("Getting sequence ...")

longN = "N"*2000000

def frowny_case(instring):
    return instring[0].lower() + instring[1:-1].upper() + instring[-1].lower()


for item in items:
    #print("length of sequence: " + str(len(out_sequence)))
    #print(item[0])
    if len(item) == 2:
        first_ctgn = last_ctgn
        lrid, last_ctgn = item
        if args.dryrun:
            lrseq = longN
        else:
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
                    if ctg["ecc"] < len(contigs[ctg["name"]]):
                        out_tmp = out_sequence[:ctg["ecc"]-len(contigs[ctg["name"]])]
                        out_sequence = out_tmp[0:-1] + out_tmp[-1].lower()
                    status = 1
                    last_used_ctg = ctg.copy()
                    # get things
                else:
                    continue
            elif status == 1:
                #tocut = len(contigs[last_used_ctg["name"]]) - (last_used_ctg["ecc"]+1) # ecc and scc are zero based
                #print("to cut: " +str(tocut) + " " + last_used_ctg["name"])
                #if tocut > 0:
                #    out_sequence = out_sequence[:-tocut]
                if last_used_ctg["ecr"] > ctg["scr"]:
                    out_sequence = out_sequence[:ctg["scr"]-last_used_ctg["ecr"]]
                else:
                    out_sequence += (lr_seq[last_used_ctg["ecr"]+1:ctg["scr"]]).lower()
                if ctg["name"] == last_ctgn:
                    out_sequence += frowny_case(contigs[ctg["name"]])
                    break
                else:
                    out_sequence += frowny_case(contigs[ctg["name"]][ctg["scc"]:ctg["ecc"]])
                last_used_ctg = ctg
    elif len(item) == 1:
        if last_ctgn == item[0]:
            continue
        else:
            #print("Th")
        # get distance between last_ctg and item[0]
            if (last_ctgn,item[0]) not in distances:
                print("get distance between " + last_ctgn + " and " + item[0])
                sys.exit()
            else:
                out_sequence = out_sequence[:distances[(last_ctgn, item[0])]]
                out_sequence += (contigs[item[0]]).upper()
                last_ctgn = item[0]

if args.logfile:
    loghandle.close()
                
if args.dryrun:
    print("Lenght of the assembly: " + str(len(out_sequence)))
else:
    with open(args.outfile, "w") as out:
        out.write(">MHC_APD\n")
        out.write(out_sequence)

