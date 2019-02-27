from argparse import ArgumentParser
from collections import defaultdict
import sys
from Bio import SeqIO


parser = ArgumentParser()
parser.add_argument("inputfile", help="Error rate file or PAF file")
parser.add_argument("cellline", help="Name of cell line")
parser.add_argument("sequencefile", help="fastq-inputfile")
parser.add_argument("contigfile", help="contig input fasta")
#parser.
parser.add_argument("outputfile", help="output in fasta")
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

contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = str(read.seq)

readids = defaultdict(set)

lreads = {}
reverse = set()
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


for rid, lr in lreads.items():
    bw = 0
    fw = 0
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.cellline): 
            if mapping["strand"] == 1:
                bw += 1
            elif mapping["strand"] == 0:
                fw += 1
            else:
                raise ValueError("strand: " + str(mapping["strand"]))
    if bw > fw:
        for mapping in lr["maps"]:
            if mapping["contig"].endswith(args.cellline): 
                mapping["strand"] = 1 if mapping["strand"] == 0 else 0
                tmp = mapping["scr"]
                mapping["scr"] = lr["length"] - mapping["ecr"]
                mapping["ecr"] = lr["length"] - tmp
                tmp = mapping["scc"]
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - tmp
        reverse.add(rid)

#a = readids["756APD"] & readids["1919APD"]
a = readids["547APD"] & readids["855APD"]
#547APD_855APD
print(a)


def get_coords(lrid, ctg1, ctg2):
    global lreads
    m1 = None
    m2 = None
    for mapping in lreads[lrid]["maps"]:
        if mapping["contig"] == ctg1:
            m1 = mapping
        if mapping["contig"] == ctg2:
            m2 = mapping
    return (m1, m2)

min_ce = -1
max_cs = -1
for lrid in a:
    #r1, r2 = get_coords(lrid, "756APD", "1919APD")
    r1, r2 = get_coords(lrid, "855APD", "547APD")
    assert(r1 != None)
    assert(r2 != None)
    if min_ce == -1 or r1['ecc'] < min_ce:
        min_ce = r1['ecc']
    if max_cs == -1 or r2['ecc'] > max_cs:
        max_cs = r2['scc']

print(min_ce)
print(max_cs)
        
    #print(c1)
    #print(c2)



#{'contig': '1919APD', 'strand': 0, 'scr': 42065, 'ecr': 46396, 'scc': 1, 'ecc': 4318, 'lenc': 4318}
#{'contig': '756APD', 'strand': 0, 'scr': 23, 'ecr': 1867, 'scc': 1027, 'ecc': 2887, 'lenc': 2902}
#{'contig': '1919APD', 'strand': 0, 'scr': 1935, 'ecr': 6311, 'scc': 5, 'ecc': 4317, 'lenc': 4318}
    

#####################
#sys.exit(0)
#####################

print("Loading sequences ... ") 

complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
def revcomp(instring):
    outstring = ""
    for char in instring:
        outstring += complement[char]
    return outstring[::-1]

lrs = {}
with open(args.sequencefile) as f:
    for i, line in enumerate(f):
        if i % 4 == 0:
            lrid = line.split(" ")[0][1:]
        elif i % 4 == 1:
            if lrid in reverse:
                lrs[lrid] = revcomp(line.rstrip())
            else:
                lrs[lrid] = line.rstrip()
print("done")


with open(args.outputfile, "w+") as out:
    for lrid in a:
        #c1, c2 = get_coords(lrid, "756APD", "1919APD")
        c1, c2 = get_coords(lrid, "855APD", "547APD")
        assert(c1 != None)
        assert(c2 != None)
        if c1['ecc'] > min_ce:
            c1['cut_left'] = c1['ecc'] - min_ce # TODO improve with Alex align tool
        else:
            c1['cut_left'] = 0

        if c2['scc'] < max_cs:
            c2['cut_right'] = max_cs - c2['scc'] # TODO improve with Alex align tool
        else:
            c2['cut_right'] = 0
        leftc = c1['ecr'] + min_ce - c1['ecc'] - 20
        rightc = c2['scr'] - max_cs + c2['scc'] + 20
        out.write(">" + lrid + "\n")
        out.write(lrs[lrid][leftc:rightc+1] + "\n")


    print(">756APD")
    print(contigs["756APD"][-40:])
    print(">1919APD")
    print(contigs["1919APD"][:40])
    print(">855APD")
    print(contigs["855APD"][-40:])
    print(">547APD")
    print(contigs["547APD"][:40])

