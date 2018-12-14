from argparse import ArgumentParser
from collections import defaultdict


parser = ArgumentParser()
parser.add_argument("inputfile", help="Error rate file or PAF file")
parser.add_argument("cellline", help="Name of cell line")
parser.add_argument("sequencefile", help="fastq-inputfile")
parser.add_argument("--paf", help="Input is paf file", action="store_true", default = False)
parser.add_argument("--blacklist", help="File containing long read ids where certain contig mappings should be ignored.")
args = parser.parse_args()



#blacklist of long reads
blacklist = {}
blacklist_contigs = set()
if args.blacklist:
    with open(args.blacklist) as f:
        for line in f:
            if line.split()[0] == "contig":
                blacklist_contigs.add(line.split()[1].rstrip())
            else:
                idx, ctg =  line.strip().split()[0:2]
                blacklist[idx] = ctg


readids = defaultdict(set)

lreads = {}
with open(args.inputfile) as f:
    for line in f:
        if args.paf:
            [rid, lenr, scr, ecr, strandstring, ctg, lenc, scc, ecc, nr_matches, block_len, quality] = line.split()[0:12]
            #print(strandstring)
            if strandstring == "+":
                strand = 0
            else:
                strand = 1
        else:
            [rid, ctg, t2, t3, t4, scr, ecr, lenr, strand, scc, ecc, lenc, t12, t13, t14, t15, t16] = line.split()
        data = {"contig":ctg,"strand":int(strand),"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
        if args.blacklist:
            if rid in blacklist:
                if blacklist[rid] == "all":
                    continue
                elif blacklist[rid] == ctg:
                    continue
            if ctg in blacklist_contigs:
                continue
        if ctg.endswith(args.cellline):
            readids[ctg].add(rid)
        if rid in lreads:
            lreads[rid]["maps"].append(data)
            if int(ecr) > lreads[rid]["rm_ecr"]:
                lreads[rid]["rm_ecr"] = int(ecr)
            if int(scr) < lreads[rid]["lm_scr"]:
                lreads[rid]["lm_scr"] = int(scr)
        else:
            lreads[rid] = {}
            lreads[rid]["length"] = int(lenr)
            lreads[rid]["maps"] = [data]
            lreads[rid]["rm_ecr"] = int(ecr)
            lreads[rid]["lm_scr"] = int(scr)

a = readids["756APD"] & readids["1919APD"]
print(a)

def revcomp(instring):
    outstring = ""
    for char in instring:
        outstring += complement[char]
    return outstring[::-1]

print("Loading sequences ... ") 
lrs = {}
with open(args.sequencefile) as f:
    for i, line in enumerate(f):
        if i % 4 == 0:
            lrid = line.split(" ")[0][1:]
        elif i % 4 == 1:
            if lrid in dot.nodes:
                if dot.nodes[lrid]["reverse"]:
                    lrs[lrid] = revcomp(line.rstrip())
                else:
                    lrs[lrid] = line.rstrip()
            elif lrid in not_merged:
                if not_merged[lrid]:
                    lrs[lrid] = revcomp(line.rstrip())
                else:
                    lrs[lrid] = line.rstrip()
    
print("done")
