from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from itertools import combinations


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")

args = parser.parse_args()

reads = {}
greads = {}
cgreads = []
contigs = {}


for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

print("Nr. of scaffolds: " + str(len(contigs)))


# nanopore read info
with open(args.efile) as f:
    for line in f:
        sline = line.split()
        rid = sline[0]
        ctg = sline[1]
        strand = int(sline[8])
        scr = int(sline[5])
        ecr = int(sline[6])
        scc = int(sline[9])
        ecc = int(sline[10])
        lc = int(sline[11])
        payload = {"contig":ctg,"strand":strand,"scr":scr,"ecr":ecr,"scc":scc,"ecc":ecc,"lc":lc}
        if rid in reads:
            reads[rid]["overlaps"].append(payload)
        else:
            reads[rid] = {}
            reads[rid]["length"] = int(sline[7])
            reads[rid]["overlaps"] = [payload]

# get interesting reads
# and sort contigs by left coordinate
greadst = {}
for rid in reads:
    counter = 0
    for item in reads[rid]["overlaps"]:
        if item["contig"].endswith("QBL"):
            greadst[rid] = reads[rid]
            break
for rid in greadst:
    #print(reads[rid]["overlaps"])
    soverlaps = sorted(reads[rid]["overlaps"], key = itemgetter("scr"))
    greads[rid] = greadst[rid]
    greads[rid]["overlaps"]=soverlaps




# cluster np-reads 
print("scaffolding long reads ....")
contig2cluster = {}
creads = {}
clusternr = 0
while len(greads) > 0:
    clusternr += 1
    current_cluster = {}
    current_contigs = set()
    # take a random read and build a cluster from it
    cr = greads.popitem()
    current_cluster[cr[0]] = cr[1]
    #print(len(current_cluster))
    olen = 0
    while len(current_cluster) != olen:
        olen = len(current_cluster)
        for contig in cr[1]["overlaps"]:
            if not contig["contig"].startswith("chr"):
                contigs.pop(contig["contig"], None)
                current_contigs.add(contig["contig"])
                contig2cluster[contig["contig"]] = clusternr
        #print("contigs: " + str(current_contigs))
        for readid,readval in greads.items():
            contig_found = False
            for contig in readval["overlaps"]:
                if contig["contig"] in current_contigs:
                    contig_found = True
                    current_cluster[readid] = readval
                    cr = (readid, greads.pop(readid))
                    break
                   
            if contig_found:
                break
    #print("cluster length: " + str(len(current_cluster)))
    creads[clusternr] = current_cluster
print("Nr. of scaffolds: " + str(clusternr+len(contigs)) + " (" + str(clusternr) + " clustered + " + str(len(contigs))+ " unclustered)")
#print(len(contig2cluster))
#print(contig2cluster)

#def lneighbour(contig):
#    print(creads[contig2cluster[contig]])

#lneighbour("2053QBL")
scaffolds={}
for i, cluster in creads.items():
    current_contigs = set([])
    for readid,read in cluster.items():
        print(read)
        for contig in read["overlaps"]:
            current_contigs.add(contig["contig"])
    scaffolds[i] = current_contigs

print(scaffolds)


sdistances = {}
with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        ctg1 = sline[0].split("_")[0].strip("+").strip("-")
        ctg2 = sline[0].split("_")[1].strip("+").strip("-")
        if sline[1] == "NA":
            continue
        moddist = float(sline[1])
        if int(ctg1.rstrip("QBL")) < int(ctg2.rstrip("QBL")):
            cstr = ctg1+"_"+ctg2
        else:
            cstr = ctg2+"_"+ctg1
        if cstr in sdistances:
            if abs(moddist) < abs(sdistances[cstr]):
                sdistances[cstr] = moddist
        else:    
            sdistances[cstr] = moddist
            
        
df = pd.DataFrame.from_dict([ldistances, sdistances])
#df.rename(index=
dc = df.T.rename(columns={0:'longread',1:'shortread'})
dc.longread = dc.longread.apply(np.mean)
dd = dc.dropna()
print(dd)

#de = dd[dd['longread']-dd['shortread']) < 30]


#get interesting differences
#print(dd[abs(dd['longread'] - dd['shortread']) > 150])


plt.scatter(dd['longread'], dd['shortread'],s= 6, alpha = 0.3)
plt.xlabel("Long Read Distances (mean: " + "{:.3f}".format(np.mean(dd['longread'])) + ")")
plt.ylabel("Short Read Distances (mean: " + "{:.3f}".format(np.mean(dd['shortread'])) + ")")

#plot line to show which values are taken
x = np.arange(-1000,1000,0.1)
y1 = x+30
y2 = x-30
plt.plot(x,y1)
plt.plot(x,y2)

plt.savefig('distances_scatter.pdf')
