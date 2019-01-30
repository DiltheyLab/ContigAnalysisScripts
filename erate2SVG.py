from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys
from itertools import combinations
from random import sample
from svgwrite.container import Group
from Bio import SeqIO



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
        data = {"contig":ctg,"strand":int(strand),"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
        if args.blacklist:
            if rid in blacklist:
                if blacklist[rid] == "all":
                    continue
                elif blacklist[rid] == ctg:
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
        if item["contig"].endswith(args.linename):
            counter += 1
    if counter >= int(args.mincontigs):
        greadst[rid] = read
        for item in read["maps"]:
            if item["contig"].endswith(args.linename):
                greads_contigset.add(item["contig"])
    #if counter == 1:
    #    intreads[rid] = read
    #    for item in read["maps"]:
    #        if item["contig"].endswith(args.linename):
    #            intreads_contigset.add(item["contig"])
            
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
                if not args.paf:
                    tmp = mapping["scc"]
                    mapping["scc"] = mapping["lenc"] - mapping["ecc"]
                    mapping["ecc"] = mapping["lenc"] - tmp
        
    # turn around and redefine wrong contigs
    for mapping in lr["maps"]:
        if mapping["contig"].endswith(args.linename): 
            if mapping["strand"] == 1: #define a new contigname and turn it around
                mapping["contig"] = mapping["contig"] + "rc"
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
whitelist = set()
clusternr = 0
if(args.whitelist):
    with open(args.whitelist) as f:
        for line in f:
            whitelist.add(line.strip())
            
    creads[0] = {}
    for rid,read in greads.items():
        for mapping in read["maps"]:
            if mapping["contig"] in whitelist:
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

csize = 0
if args.alignreads:
    clustered_reads = sorted_reads
else:
    clustered_reads = creads
clip_path_idc = 0
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
            ctg = read["contig"].rstrip(args.linename)
            ctgn = read["contig"].rstrip("rc")
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctg = ctg[0:8]
                continue

            clip_path_idc+=1
            clip_path = dwg.defs.add(dwg.clipPath(id=str(clip_path_idc)))
            clip_path.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100), ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1))
            leftclip = scc/100
            rightclip = (contigs[ctgn]-ecc)/100
            #g.add(svgwrite.shapes.Rect((xpad+(sc/100)-leftflip,yoff+ypos-ctg_y_halfdrawsize), ((ec-sc)/100+leftclip+rightclip,ctg_y_drawsize), stroke='black', stroke_width=1, fill = 'url(#rwb)'))
            
            g.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100)-leftclip,ypad+ypos-6), ((ec-sc)/100+leftclip+rightclip,12), stroke='black', stroke_width=1, fill='url(#rwb)', clip_path='url(#'+str(clip_path_idc)+')'))
            if above:
                yt = ypad+ypos-8
                if col == "blue":
                    col = "black"
                else: 
                    col = "blue"
            else:
                yt = ypad+ypos+13
            above = not above
            #g.add(dwg.text("s"+ctg+"$", insert=(xpad+(xoffset+sc)/100,yt),fill=col, style="font-size:4"))
            g.add(dwg.text(ctg, insert=(xpad+(xoffset+sc)/100,yt),fill=col, style="font-size:6"))
            if read["strand"] == 0:
                direction = ">"
            elif read["strand"] == 1:
                direction = "<"
            g.add(dwg.text(direction, insert=(xpad+(xoffset+sc)/100,ypad+ypos+2),style="font-size:6"))
        dwg.add(g)
    ypos += 28
    dwg.add(dwg.line((10, ypad+ypos), ( 10000, ypad+ypos+10), stroke=svgwrite.rgb(0, 0, 0, '%')))


dwg.save()


