from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys
from itertools import combinations
from random import sample
from svgwrite.container import Group
from Bio import SeqIO
from collections import defaultdict

from scaffold import Scaffold


parser = ArgumentParser()
parser.add_argument("inputfile", help="Error-rate or paf-file")
parser.add_argument("contigfile", help="Fasta File containing contigs")
parser.add_argument("--paf", help="Input is paf file", action="store_true", default = False)
parser.add_argument("--whitelist", help="Only plot long reads with ids found in this whitelist file.")
parser.add_argument("--blacklist", help="Do not plot long reads in this file.")
parser.add_argument("--alignreads", action="store_true", help="Reads will be turned around if necessary and aligned according to their distances.")
parser.add_argument("--mincontigs", default=2, help="Minimum number of contigs that have to map into the long read, for the long read to be considered.")
parser.add_argument("linename", help="Name of the cellline")
parser.add_argument("output", help="SVG output file")

args = parser.parse_args()

lreads ={}
greads = {}
cgreads = []

#blacklist of long reads
blacklist = defaultdict(set)
blacklist_contigs = set()
if args.blacklist:
    with open(args.blacklist) as f:
        for line in f:
            if line.split()[0] == "contig":
                blacklist_contigs.add(line.split()[1].rstrip())
            else:
                idx, ctg =  line.strip().split()[0:2]
                blacklist[idx].add(ctg)


contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

# paf format 
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
        data = {"name":ctg,"strand":int(strand),"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
        if args.blacklist:
            if rid in blacklist:
                if "all" in blacklist[rid]:
                    continue
                elif ctg in blacklist[rid]:
                    continue
            if ctg in blacklist_contigs:
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
        if item["name"].endswith(args.linename):
            counter += 1
    if counter >= int(args.mincontigs):
        greadst[rid] = read
        for item in read["maps"]:
            if item["name"].endswith(args.linename):
                greads_contigset.add(item["name"])
    #if counter == 1:
    #    intreads[rid] = read
    #    for item in read["maps"]:
    #        if item["name"].endswith(args.linename):
    #            intreads_contigset.add(item["name"])
            
#intersection = greads_contigset.intersection(intreads_contigset)
#sintersection = sorted(intersection, key = lambda x: int(x.rstrip(args.linename)))
#print(sintersection)


for rid in intreads:
    if len(intreads[rid]["maps"]) == 1:
        pass
        #print("\t".join([rid, str(intreads[rid])]))
        

greadst.update(intreads)

# turn reads around if necessary
for rid, lr in greadst.items():
    bw = 0
    fw = 0
    for mapping in lr["maps"]:
        if mapping["name"].endswith(args.linename): 
            if mapping["strand"] == 1:
                bw += 1
            elif mapping["strand"] == 0:
                fw += 1
            else:
                raise ValueError("strand: " + str(mapping["strand"]))
    if bw > fw:
        for mapping in lr["maps"]:
            if mapping["name"].endswith(args.linename): 
                mapping["strand"] = 1 if mapping["strand"] == 0 else 0
                tmp = mapping["scr"]
                mapping["scr"] = lr["length"] - mapping["ecr"]
                mapping["ecr"] = lr["length"] - tmp
                if not args.paf:
                    tmp = mapping["scc"]
                    mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                    mapping["ecc"] = mapping["lenc"] - tmp
        
    # turn around and redefine wrong contigs
    for mapping in lr["maps"]:
        if mapping["name"].endswith(args.linename): 
            if mapping["strand"] == 1: #define a new contigname and turn it around
                mapping["name"] = mapping["name"] + "rc"
                if not args.paf:
                    tmp = mapping["scc"]
                    mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                    mapping["ecc"] = mapping["lenc"] - tmp
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
ccontigs = {}
whitelist = set()
clusternr = 0
if(args.whitelist):
    with open(args.whitelist) as f:
        for line in f:
            whitelist.add(line.strip())
            
    creads[0] = {}
    for rid,read in greads.items():
        for mapping in read["maps"]:
            if mapping["name"] in whitelist:
                creads[0][rid] = read
else:
    while len(greads) > 0:
        clusternr += 1
        current_cluster = {}
        current_contigs = set()
        # take a random read and build a cluster from it

        ck = sample(greads.keys(),1)[0]
        #print(ck)
        cr = (ck, greads[ck])
        #cr = greads.popitem()
        current_cluster[cr[0]] = cr[1]
        #for contig in cr[1]["maps"]:
        #    if not contig["name"].startswith("chr"):
        #        current_contigs.add(contig["name"])
        #print(len(current_cluster))
        olen = 0
        while len(current_cluster) != olen:
            olen = len(current_cluster)
            for contig in cr[1]["maps"]:
                if not contig["name"].startswith("chr"):
                    current_contigs.add(contig["name"])
            del(greads[cr[0]])
            #print("contigs: " + str(current_contigs))
            for readid,readval in greads.items():
                contig_found = False
                for contig in readval["maps"]:
                    if contig["name"] in current_contigs:
                        contig_found = True
                        current_cluster[readid] = readval
                        cr = (readid, greads[readid])
                        break
                if contig_found:
                    break
            
        print("cluster length: " + str(len(current_cluster)))
        creads[clusternr] = current_cluster
        ccontigs[clusternr] = current_contigs

def compare_longreads(lr1, lr2):
    l1c = []
    l2c = []
    for m in lr1["maps"]:
        cn = m["name"]
        if not cn.startswith("chr"):
            l1c.append(cn)
    for m in lr2["maps"]:
        cn = m["name"]
        if not cn.startswith("chr"):
            l2c.append(cn)
    common_ctgs = set(l1c).intersection(set(l2c))
    return common_ctgs


def get_contig_info(lr,ctg):
    for maps in lr["maps"]:    
        if maps["name"] == ctg:
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

if args.alignreads:
    for lrs in combinations(ogreads.keys(), 2):
        lr1 = ogreads[lrs[0]]
        lr2 = ogreads[lrs[1]]
        common_ctgs = compare_longreads(lr1, lr2)
        if len(common_ctgs) > 0:
            dists = get_distances(lr1, lr2, common_ctgs)
            lr_dists[lrs[0]][lrs[1]]=dists[0]
            lr_dists[lrs[1]][lrs[0]]=-dists[0]

                
sorted_reads = {}
if args.alignreads:
    for cluster in creads:
        #print("ccontigs: " + str(ccontigs[cluster]))
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
        #print("A: "+ str(A))
        #print("B: "+ str(B))
        U = A-B
        #print("U: "+ str(U))
        UC = U.copy()
        while True:
        #for i in range(100):
            #print(U)
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
        
        sorted_reads[cluster] = sorted(creads[cluster].keys(), key=lambda x: lr_dists[srid][x])


        # store resulsts of this
        creads[cluster]["smallest_id"] = srid #for this id all values are present in lr_dists
        creads[cluster]["smallest_offset"] = -smval

#sys.exit(0)


# draw interesting reads
ypos = 0
xtext = 10
xpad = 200
ypad = 10

dwg = svgwrite.Drawing(args.output, size=(u'4000', u'10000'), profile='full')
lineargrad = dwg.defs.add(svgwrite.gradients.LinearGradient(id="rwb"))
lineargrad.add_stop_color("0%","#FF0000")
lineargrad.add_stop_color("50%","#FFFFFF")
lineargrad.add_stop_color("100%","#0000FF")


def shortname(ctgname):
    if "_" in ctgname:
        return ctgname.split("_")[1]
    else:
        return ctgname

csize = 0
if args.alignreads:
    clustered_reads = sorted_reads
else:
    clustered_reads = creads
gradient_idc = 0
for cluster in clustered_reads:
    for rid in clustered_reads[cluster]:
        if rid.startswith("small"):
            continue
        if args.alignreads:
            xoffset = creads[cluster]["smallest_offset"] + lr_dists[creads[cluster]["smallest_id"]][rid]
        else:
            xoffset = 0
        ypos += 28
        dwg.add(dwg.text(rid, insert=(xtext, ypad+ypos+1), fill='black', style="font-size:7"))
        g = dwg.defs.add(Group(id=rid))
        g.add(dwg.line((xpad+ xoffset/100, ypad+ypos), ( xpad + (xoffset+ogreads[rid]["length"])/100, ypad+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))
        above = True
        col = "black"
        for read in ogreads[rid]["maps"]:
            #print(read)
            sc = read["scr"]
            ec = read["ecr"]
            scc = read["scc"]
            ecc = read["ecc"]
            ctg = read["name"].rstrip("rc").rstrip(args.linename)
            ctgn = read["name"].rstrip("rc")
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctg = ctg[0:8]
                continue
            lenc = contigs[shortname(ctgn)]

            gradient_idc += 1
            #clip_path = dwg.defs.add(dwg.clipPath(id=str(clip_path_idc)))
            #clip_path.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100), ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1))
            #leftclip = scc/100
            #rightclip = (contigs[ctgn]-ecc)/100
            lineargrad = dwg.defs.add(svgwrite.gradients.LinearGradient(id=str(gradient_idc), x1=-scc/(ecc-scc), x2=1+(contigs[shortname(ctgn)]-ecc)/(ecc-scc), y1=0, y2=0))
            if ctg.endswith("0") or ctg.endswith("1"):
                col1 = "#FF0000"
                col2 = "#0000FF"
            elif ctg.endswith("2") or ctg.endswith("3"):
                col1 = "#FFE119"
                col2 = "#000000"
            elif ctg.endswith("4") or ctg.endswith("5"):
                col1 = "#911eb4"
                col2 = "#a9a9a9"
            elif ctg.endswith("6") or ctg.endswith("7"):
                col1 = "#000075"
                col2 = "#f58231"
            elif ctg.endswith("8") or ctg.endswith("9"):
                col1 = "#e6beff"
                col2 = "#808000"
                
            lineargrad.add_stop_color("0%",col1)
            lineargrad.add_stop_color("50%","#FFFFFF")
            lineargrad.add_stop_color("100%",col2)
            #g.add(svgwrite.shapes.Rect((xpad+(sc/100)-leftflip,yoff+ypos-ctg_y_halfdrawsize), ((ec-sc)/100+leftclip+rightclip,ctg_y_drawsize), stroke='black', stroke_width=1, fill = 'url(#rwb)'))
#lineargrad = dwg.defs.add(svgwrite.gradients.LinearGradient(id="rwb"))
            g.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100),ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1, fill='url(#'+str(gradient_idc)+')'))
            if above:
                yt = ypad+ypos-8
                if col == "blue":
                    col = "black"
                else: 
                    col = "blue"
            else:
                yt = ypad+ypos+13
            above = not above
            #g.add(dwg.text(ctg, insert=(xpad+(xoffset+sc)/100,yt),fill=col, style="font-size:6"))
            g.add(dwg.text("s"+ctg+"$", insert=(xpad+(xoffset+sc)/100,yt),fill=col, style="font-size:4"))
            if read["strand"] == 0:
                direction = ">"
            elif read["strand"] == 1:
                direction = "<"
            g.add(dwg.text(direction, insert=(xpad+(xoffset+sc)/100,ypad+ypos+2),style="font-size:6"))
        dwg.add(g)
    ypos += 28
    dwg.add(dwg.line((10, ypad+ypos), ( 10000, ypad+ypos+10), stroke=svgwrite.rgb(0, 0, 0, '%')))


dwg.save()


