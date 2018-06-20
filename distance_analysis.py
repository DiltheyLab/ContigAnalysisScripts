from argparse import ArgumentParser
from operator import itemgetter
import sys


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="contig distance summary file")

args = parser.parse_args()


reads = {}
greads = {}
cgreads = []

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
greadst = {}
for rid in reads:
    counter = 0
    for item in reads[rid]["overlaps"]:
        if item["contig"].endswith("QBL"):
            counter += 1
    if counter >= 2:
        greadst[rid] = reads[rid]
            
# sort contigs by left coordinate
for rid in greadst:
    #print(reads[rid]["overlaps"])
    soverlaps = sorted(reads[rid]["overlaps"], key = itemgetter("scr"))
    greads[rid] = greadst[rid]
    greads[rid]["overlaps"]=soverlaps

#print(greads)

distances = {}
for rid in greads:
    oviter = iter(greads[rid]["overlaps"])
    #print(ov)
    try:
        ovold = next(oviter)
        if ovold["contig"].startswith("chr"):
            continue
    except StopIteration:
        continue
    while True:
        try:
            ovnew = next(oviter)
            if ovnew["contig"].startswith("chr"):
                continue
        except StopIteration:
            break
        #print("distance between " + ovold["contig"] + " and " + ovnew["contig"] + ": " + str(ovnew["scr"] - ovold["ecr"]- ovold["lc"] + ovold["ecc"] - ovnew["scc"]))
        if ovnew["contig"] == ovold["contig"]:
            continue
        if ovnew["strand"] == 0 and ovold["strand"] == 0:
            #cstring = ovold["contig"] + "_" + ovnew["contig"]
            distance = ovnew["scr"] - ovold["ecr"]- (ovold["lc"] - ovold["ecc"]) - ovnew["scc"] + 1
        elif ovnew["strand"] == 1 and ovold["strand"] == 1:
            #cstring = ovnew["contig"] + "_"  + ovold["contig"]
            distance = ovnew["scr"] - ovold["ecr"]- (ovold["lc"] - ovold["ecc"]) - ovnew["scc"] + 1
        else:
            #print(str(ovold)+ " " + str(ovnew))
            continue
        if int(ovold["contig"].rstrip("QBL")) < int(ovnew["contig"].rstrip("QBL")):
            cstring = ovold["contig"] + "_" + ovnew["contig"]
        else:
            cstring = ovnew["contig"] + "_" + ovold["contig"]
    
        if cstring in distances:
            distances[cstring].append(distance)
        else:
            distances[cstring] = [distance]
        ovold = ovnew

        
for key, value in distances.items():
    #print(str(key) + " " + str(value)) 
    if len(value)>1:
        print(str(key) + "\t" + str(value)) 

distances2 = {}

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
        if cstr in distances2:
            if abs(moddist) < abs(distances2[cstr]):
                distances2[cstr] = moddist
        else:    
            distances2[cstr] = moddist
            
for pair in distances2:
    if pair in distances:
        npdists = distances[pair]
        #print(pair + "\t" + str(sum(npdists)/float(len(npdists))) + "\t" + str(distances2[pair]))
        
