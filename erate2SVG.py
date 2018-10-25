from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys
from itertools import combinations
from random import sample



parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("--whitelist", help="Only plot long reads with ids found in this whitelist file.")
parser.add_argument("--blacklist", help="Do not plot long reads in this file.")
parser.add_argument("--alignreads", action="store_true", help="Reads will be turned around if necessary and aligned according to their distances.")
parser.add_argument("linename", help="Name of the cellline")
parser.add_argument("output", help="SVG output file")

args = parser.parse_args()

lreads ={}
greads = {}
cgreads = []

#blacklist of long reads
blacklist = {}
if args.blacklist:
    with open(args.blacklist) as f:
        for line in f:
            idx, ctg =  line.strip().split()[0:2]
            blacklist[idx] = ctg

# nanopore read info
with open(args.efile) as f:
    for line in f:
        [rid, ctg, t2, t3, t4, scr, ecr, lenr, strand, scc, ecc, lenc, t12, t13, t14, t15, t16] = line.split()
        data = {"contig":ctg,"strand":int(strand),"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
        if args.blacklist:
            if rid in blacklist:
                if blacklist[rid] == "all":
                    continue
                elif blacklist[rid] == ctg:
                    continue
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

# get interesting reads
greadst = {}
intreads = {}
greads_contigset = set()
intreads_contigset = set()

for rid, read in lreads.items():
    counter = 0
    for item in read["maps"]:
        if item["contig"].endswith(args.linename):
            counter += 1
    if counter >= 2:
        greadst[rid] = read
        for item in read["maps"]:
            if item["contig"].endswith(args.linename):
                greads_contigset.add(item["contig"])
    if counter == 1:
        intreads[rid] = read
        for item in read["maps"]:
            if item["contig"].endswith(args.linename):
                intreads_contigset.add(item["contig"])
            
intersection = greads_contigset.intersection(intreads_contigset)
sintersection = sorted(intersection, key = lambda x: int(x.rstrip(args.linename)))
#print(sintersection)


for rid in intreads:
    if len(intreads[rid]["maps"]) == 1:
        pass
        #print("\t".join([rid, str(intreads[rid])]))
        
# turn reads around

# turn reads around if necessary
for rid, lr in greadst.items():
    bw = 0
    fw = 0
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.linename): 
            if mapping["strand"] == 1:
                bw += 1
            elif mapping["strand"] == 0:
                fw += 1
            else:
                raise ValueError("strand: " + str(mapping["strand"]))
    if bw > fw:
        for mapping in lr["maps"]:
            if mapping["contig"].endswith(args.linename): 
                mapping["strand"] = 1 if mapping["strand"] == 0 else 0
                tmp = mapping["scr"]
                mapping["scr"] = lr["length"] - mapping["ecr"]
                mapping["ecr"] = lr["length"] - tmp
                tmp = mapping["scc"]
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - tmp
        
    # turn around and redefine wrong contigs
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.linename): 
            if mapping["strand"] == 1: #define a new contigname and turn it around
                mapping["contig"] = mapping["contig"] + "rc"
                mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                mapping["ecc"] = mapping["lenc"] - mapping["scc"]
                mapping["strand"] = 0


# sort contigs by left coordinate
for rid in greadst:
    soverlaps = sorted(lreads[rid]["maps"], key = itemgetter("scr"))
    greads[rid] = greadst[rid]
    greads[rid]["maps"]=soverlaps

ogreads = greads.copy()
print("length greads: " + str(len(greads)))


# cluster np-reads or only keep whitelisted ones
creads = {}
whitelist = set()
clusternr = 0
if(args.whitelist):
    with open(args.whitelist) as f:
        for line in f:
            whitelist.add(line.strip())
            
    creads[0] = {}
    for rid,read in greads.items():
        for mapping in read["maps"]:
            if mapping[0] in whitelist:
                creads[0][rid] = read
else:
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
            for contig in cr[1]["maps"]:
                if not contig["contig"].startswith("chr"):
                    current_contigs.add(contig["contig"])
            #print("contigs: " + str(current_contigs))
            for readid,readval in greads.items():
                contig_found = False
                for contig in readval["maps"]:
                    if contig["contig"] in current_contigs:
                        contig_found = True
                        current_cluster[readid] = readval
                        cr = (readid, greads.pop(readid))
                        break
                if contig_found:
                    break
            
        print("cluster length: " + str(len(current_cluster)))
        creads[clusternr] = current_cluster

def compare_longreads(lr1, lr2):
    l1c = []
    l2c = []
    for m in lr1["maps"]:
        cn = m["contig"]
        if cn.endswith(args.linename):
            l1c.append(cn)
    for m in lr2["maps"]:
        cn = m["contig"]
        if cn.endswith(args.linename):
            l2c.append(cn)
    common_ctgs = set(l1c).intersection(set(l2c))
    return common_ctgs


def get_contig_info(lr,ctg):
    for maps in lr["maps"]:    
        if maps["contig"] == ctg:
            return maps

# calculate distances for each cluster
def get_distances(lr1,lr2, common_ctgs):
    dists = []
    for ctg in common_ctgs:
        m1 = get_contig_info(lr1,ctg)
        m2 = get_contig_info(lr2,ctg)
        dists.append((m1["scr"]-m1["scc"]) - (m2["scr"]-m2["scc"]))
    return dists
            
#lr_dist = defaultdict(lambda:([],[]))
lr_dists = {}
#initialize matrix
for lid, read in lreads.items():
    lr_dists[lid] = {lid:0}

for lrs in combinations(ogreads.keys(), 2):
    lr1 = ogreads[lrs[0]]
    lr2 = ogreads[lrs[1]]
    common_ctgs = compare_longreads(lr1, lr2)
    if len(common_ctgs) > 0:
        dists = get_distances(lr1, lr2, common_ctgs)
        lr_dists[lrs[0]][lrs[1]]=dists[0]
        lr_dists[lrs[1]][lrs[0]]=-dists[0]

            
for cluster in creads:
    #print(creads[cluster].keys())
    # find leftmost read 
    while True:
        rid = sample(creads[cluster].keys(),1)[0]
        for rid2 in creads[cluster].keys():
            if rid2 in lr_dists[rid]:
                if lr_dists[rid][rid2] < 0:
                    rid = rid2
                    break
        break
    srid = rid

    # fill lr_dists with transitive distances
    A = set(creads[cluster].keys())
    B = set(lr_dists[srid].keys())
    U = A-B
    UC = U.copy()
    while True:
        for rid in U:
            for rid2 in lr_dists[rid]:
                if rid2 not in UC:
                    lr_dists[srid][rid] = lr_dists[srid][rid2] + lr_dists[rid2][rid]
                    lr_dists[rid][srid] = lr_dists[srid][rid]
                    UC.remove(rid)
                    break
        U = UC.copy()
        if not U:
            break
   
    # now there is a new leftmost read
    smval = 0
    nsrid = srid
    for rid2 in lr_dists[srid]:
        if lr_dists[srid][rid2] < smval:
            smval = lr_dists[srid][rid2]
            nsrid = rid2

    # rid should be smallest
    creads[cluster]["smallest_id"] = srid
    creads[cluster]["smallest_offset"] = -smval
#sys.exit(0)


# draw interesting reads
ypos = 0
xtext = 10
xpad = 200
ypad = 10
dwg = svgwrite.Drawing(args.output, size=(u'1700', u'4600'), profile='full')
csize = 0
for cluster in creads:
    for rid,read in creads[cluster].items():
        if rid.startswith("small"):
            continue
        if args.alignreads:
            xoffset = creads[cluster]["smallest_offset"] + lr_dists[creads[cluster]["smallest_id"]][rid]
        else:
            xoffset = 0
        ypos += 20
        dwg.add(dwg.text(rid, insert=(xtext, ypad+ypos+1), fill='black', style="font-size:7"))
        dwg.add(dwg.line((xpad+ xoffset/100, ypad+ypos), ( xpad + (xoffset+ogreads[rid]["length"])/100, ypad+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))
        above = True
        col = "black"
        for read in ogreads[rid]["maps"]:
            #print(read)
            sc = read["scr"]
            ec = read["ecr"]
            ctg = read["contig"].rstrip(args.linename)
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctg = ctg[0:8]
            dwg.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100),ypad+ypos-3), ((ec-sc)/100,6), stroke='black', stroke_width=1, fill = 'white'))
            if above:
                yt = ypad+ypos-4
                if col == "blue":
                    col = "black"
                else: 
                    col = "blue"
            else:
                yt = ypad+ypos+7
            above = not above
            dwg.add(dwg.text(ctg, insert=(xpad+(xoffset+sc)/100,yt),fill=col, style="font-size:4"))
            if read["strand"] == 0:
                direction = ">"
            elif read["strand"] == 1:
                direction = "<"
            dwg.add(dwg.text(direction, insert=(xpad+(xoffset+sc)/100,ypad+ypos+2),style="font-size:6"))
    dwg.add(dwg.line((10, ypad+ypos+10), ( 1000, ypad+ypos+10), stroke=svgwrite.rgb(0, 0, 0, '%')))


dwg.save()





