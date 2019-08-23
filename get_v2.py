from argparse import ArgumentParser
from Bio import SeqIO
import pysam
import sys
import re
from collections import defaultdict
from scaffold import Longreads, LongReadSVG, revcomp

parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("sequencefile", help="Input File in BAM format")
parser.add_argument("contigstringfile", help="Input File containing subsequent longread ids")
parser.add_argument("contigfile", help="Contig File in FASTA format")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("--logfile", help="Logging File")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("outfile", help="Sequence Output File")
args = parser.parse_args()

# Get Contig Strings
contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = str(read.seq)

# Get Blacklist
blacklist = defaultdict(list)
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist[sline[1]] = float(sline[2])
            else:
                blacklist[sline[0]].append(sline[1])

def get_query_pos(record, ref_position):
    rpos = record.reference_start
    qpos = 0
    for opt in record.cigartuples:
        #print("rpos: " + str(rpos))
        op, times = opt
        #print(opt)
        if op in [pysam.CINS, pysam.CSOFT_CLIP]: # consumes query
            qpos += times
        elif op == pysam.CHARD_CLIP:
            continue
        elif op == pysam.CDEL:
            if rpos + times > ref_position:
                break
            rpos += times
        elif op == pysam.CMATCH:
            if rpos + times > ref_position:
                qpos += ref_position - rpos
                break
            else:
                rpos += times
                qpos += times
        else: 
            print("operation not implemented: " + op)
            sys.exit()
    return qpos
    

# Get Bamfile
bamfile = pysam.AlignmentFile(args.sequencefile,"rb")
#rstart = 125 
#rstop = 145 
iters = bamfile.fetch("MHC_APD_v1c",rstart,rstop)
for rec in iters:
    print(str(rec).split()[0:5])
    qstart = get_query_pos(rec, rstart)
    qstop = get_query_pos(rec, rstop)
    print("qstart:qstop [" + str(qstart) + ":" + str(qstop) + "]")
    if rec.query_sequence:
        print(rec.query_sequence[qstart:qstop])
    #print("rpos: " + str(rpos) + " qpos: " + str(qpos))
         
#sys.exit()
    

# Logfile
if args.logfile:
    loghandle = open(args.logfile,"w+")

print("Nr. of Contigs: " + str(len(contigs)))

scafs = Longreads(args.inputfiles, blacklist, args.linename)
#scafs.filter_low_quality_contigs(0.76)
scafs.turn_longreads_around(logging=loghandle)
scafs.sort_by_starts()
scafs.filter_contigs_by_coverage(0.5,ignore_ends=True, verbose = False)

print("Nr. of reads: " + str(len(scafs.lreads)))

# Parse Recipe File
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

# Check that fastq is present
for item in items:
    if len(item) == 2:
        if item[0] not in lreadseqs:
            print(str(item[0]) +" not found in sequence file. Aborting !")
            sys.exit()

# some curated distances
distances = {}
distances[("2000APD","1503APD")] = -20
distances[("1503APD","476APD")] = -13
distances[("476APD","1405APD")] = -9
distances[("2316APD","1498APD")] = -12
distances[("504APD","49APD")] = -447
distances[("849APD","1471APD")] = -11
distances[("1192DBB","462DBB")] = -9712

# Get first contig
firstcontig = items[0][0]
print("First contig: " + firstcontig)
items = items[1:]

out_sequence = contigs[firstcontig]
overall_length = len(out_sequence)
last_ctgn = firstcontig

# Check that important contigs weren't filtered out
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

print("Getting sequence ...")

longN = "N"*2000000

for item in items:
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
                #tocut = len(contigs[last_used_ctg["name"]]) - (last_used_ctg["ecc"]+1) # ecc and scc are zero based
                #print("to cut: " +str(tocut) + " " + last_used_ctg["name"])
                #if tocut > 0:
                #    out_sequence = out_sequence[:-tocut]
                if last_used_ctg["ecr"] > ctg["scr"]:
                    out_sequence = out_sequence[:ctg["scr"]-last_used_ctg["ecr"]]
                else:
                    out_sequence += lr_seq[last_used_ctg["ecr"]+1:ctg["scr"]]
                out_sequence += contigs[ctg["name"]][ctg["scc"]:ctg["ecc"]]
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

if args.logfile:
    loghandle.close()
                
else:
    with open(args.outfile, "w") as out:
        out.write(">MHC_APD\n")
        out.write(out_sequence)

