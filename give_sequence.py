from argparse import ArgumentParser
from Bio import SeqIO
from scaffold import Longreads, revcomp
from collections import defaultdict


parser = ArgumentParser()
parser.add_argument("inputfiles", help="Input Files in Error-Rate or PAF format", nargs="+")
parser.add_argument("sequencefile", help="Sequences of specific reads")
parser.add_argument("contigs", help="Specify the contigs. Format: '3APD-823APD'")
parser.add_argument("lrid", help="Specify the longread id")
parser.add_argument("linename", help="Line name")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
args = parser.parse_args()

blacklist = defaultdict(list)
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist[sline[1]] = "y"
            else:
                blacklist[sline[0]].append(sline[1])


#print(set([args.lrid]))
scafs = Longreads(args.inputfiles, blacklist, args.linename, whitelist_lreads=set([args.lrid]))
scafs.turn_longreads_around()
scafs.sort_by_starts()
#print(scafs.lreads)

lrseqs = dict()
for read in SeqIO.parse(args.sequencefile, "fastq"):
    lrseqs[read.id] = str(read.seq)

ctg1, ctg2 = args.contigs.split("-")
read = scafs.lreads[args.lrid]
sc = 0
ec = 0
for ctg in read["maps"]:
    if ctg["name"] == ctg1:
        #print("\t".join([ctg["name"], str(ctg["scr"]), str(ctg["ecr"]), str(ctg["scc"]), str(ctg["ecc"])]))
        sc = ctg["ecr"]
    if ctg["name"] == ctg2:
        ec = ctg["scr"]-1

print("start: " + str(sc))
print("end: " + str(ec))

if "reverse" in read:
    tmp = ec
    ec = read["length"] - sc
    sc = read["length"] - tmp
    print("start: " + str(sc))
    print("end: " + str(ec))
    print(revcomp(lrseqs[args.lrid][sc:ec]))
else:
    print(lrseqs[args.lrid][sc:ec])



