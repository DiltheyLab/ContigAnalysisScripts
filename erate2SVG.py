from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys
from itertools import combinations
from random import sample
from svgwrite.container import Group
from Bio import SeqIO
from collections import defaultdict, deque
from scaffold import Scaffold, Scaffolds


parser = ArgumentParser()
parser.add_argument("inputfiles", help="Error-rate or paf-file", nargs="+")
parser.add_argument("contigfile", help="Fasta File containing contigs")
#parser.add_argument("--paf", help="Input is paf file", action="store_true", default = False)
parser.add_argument("--whitelist", help="Only plot long reads with ids found in this whitelist file.")
parser.add_argument("--whitelist_lrs", help="Only plot long reads with these ids.")
parser.add_argument("--blacklist", help="Do not plot long reads in this file.")
parser.add_argument("--alignreads", action="store_true", help="Reads will be turned around if necessary and aligned according to their distances.")
parser.add_argument("--mincontigs", default=2, help="Minimum number of contigs that have to map into the long read, for the long read to be considered.")
parser.add_argument("linename", help="Name of the cellline")
parser.add_argument("output", help="SVG output file")

args = parser.parse_args()

# blacklist of long reads
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

whitelist_ctgs = set()
whitelist_lreads = set()
if(args.whitelist):
    with open(args.whitelist) as f:
        for line in f:
            whitelist_ctgs.add(line.strip())
elif(args.whitelist_lrs):
    with open(args.whitelist_lrs) as f:
        for line in f:
            whitelist_lreads.add(line.strip())

contigs = {}
for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)


scafs = Scaffolds(args.inputfiles, blacklist, args.linename, whitelist_lreads)
if whitelist_ctgs:
    scafs.filter_whitelist_ctgs(whitelist_ctgs)
scafs.filter_contigcounts(args.mincontigs)
scafs.turn_longreads_around()
scafs.sort_by_starts()


print("Reads meeting criteria: " + str(len(scafs.lreads)))
lr_dists = scafs.pseudoalign_all()


def nodes_fully_connected(nodes):
    for node in nodes:
        if nodes - lr_dists[node].keys():
            return False
    else:
        return True

def get_least_connected_node(nodes):
    least_connections = len(nodes)
    least_node = None
    for node in nodes:
        if len(lr_dists[node].keys() & nodes) < least_connections:
            least_connections = len(lr_dists[node])
            least_node = node
    return least_node

all_nodes = set(scafs.lreads.keys())
creads = []
while all_nodes:
    nodes = all_nodes.copy()
    while not nodes_fully_connected(nodes):
        nos = get_least_connected_node(nodes)
        nodes.remove(nos)
    origin = sample(nodes, 1)[0]
    creads.append(deque(sorted(nodes, key=lambda x: lr_dists[origin][x])))
    all_nodes -=nodes



#lreads = scafs.lreads.copy()


# find clusters:
#creads = []
#while lreads:
    # build cluster from random chosen read
#    lrid, lread = lreads.popitem()
#    cluster = deque([lrid])
#    to_analyze = deque([lrid])
#    while to_analyze:
#        citem = to_analyze.popleft()
#        nitems = set(lr_dists[citem].keys()) & lreads.keys()
#        for item in nitems:
#            cluster.append(item)
#            to_analyze.append(item)
#            del(lreads[item])
#    creads.append(cluster)
print("Number of clusters found: " + str(len(creads)))
for item in creads:
    print("Number of reads: " + str(len(item)))

sorted_reads = []
smids = []
smoffs = []
sorted_clusters = []
if args.alignreads:
    for cluster in creads:
        clst = cluster.copy()
        origin = clst.popleft()
        distance_known = set([origin])
        while clst:
            citem = clst.popleft()
            if citem not in lr_dists[origin]:
                for bitem in distance_known:
                    if bitem in lr_dists[citem]:
                        lr_dists[origin][citem] = lr_dists[origin][bitem] + lr_dists[bitem][citem]
                        lr_dists[citem][origin] = -(lr_dists[origin][bitem] + lr_dists[bitem][citem])
            distance_known.add(citem)
        #print(lr_dists[origin])
        #print(len(lr_dists[origin]))
        sorted_reads = sorted(set(lr_dists[origin].keys()) & set(cluster), key = lambda x: lr_dists[origin][x])
        smid = sorted_reads[0]
        #print(smid)
        for item in lr_dists[origin]:
            if item not in lr_dists[smid]:
                lr_dists[smid][item] = lr_dists[smid][origin] + lr_dists[origin][item]
        #print(lr_dists[smid])
        #print(len(lr_dists[smid]))

        sorted_clusters.append(sorted_reads)



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
#mask1 = dwg.defs.add(svgwrite.masking.Mask(maskContentUnits = 'objectBoundingBox', id="mask1"))
#mask1.add(svgwrite.shapes.Rect(fill="black", x="0", y="0", width="100%", height="100%"))
#mask1.add(svgwrite.path.Path(d="M 0 0 L 1 1 L 0 1", fill="white"))
#mask2 = dwg.defs.add(svgwrite.masking.Mask(maskContentUnits = 'objectBoundingBox', id="mask2"))
#mask2.add(svgwrite.shapes.Rect(fill="black", x="0", y="0", width="100%", height="100%"))
#mask2.add(svgwrite.path.Path(d="M 0 0 L 1 1 L 1 0", fill="white"))



def shortname(ctgname):
    if "_" in ctgname:
        return ctgname.split("_")[1]
    else:
        return ctgname

csize = 0
if not args.alignreads:
    sorted_clusters = creads

gradient_idc = 0
for cluster in sorted_clusters:
    for rid in cluster:
        if args.alignreads:
            xoffset = lr_dists[cluster[0]][rid]
        else:
            xoffset = 0
        ypos += 28
        dwg.add(dwg.text(rid, insert=(xtext, ypad+ypos+1), fill='black', style="font-size:7"))
        g = dwg.defs.add(Group(id=rid))
        g.add(dwg.line((xpad+ xoffset/100, ypad+ypos), ( xpad + (xoffset+scafs.lreads[rid]["length"])/100, ypad+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))
        above = True
        col = "black"
        for read in scafs.lreads[rid]["maps"]:
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
                
            if read["strand"] == 1:
                tcol = col1
                col1 = col2
                col2 = tcol
            lineargrad.add_stop_color("0%",col1)
            lineargrad.add_stop_color("50%","#FFFFFF")
            lineargrad.add_stop_color("100%",col2)
            #g.add(svgwrite.shapes.Rect((xpad+(sc/100)-leftflip,yoff+ypos-ctg_y_halfdrawsize), ((ec-sc)/100+leftclip+rightclip,ctg_y_drawsize), stroke='black', stroke_width=1, fill = 'url(#rwb)'))
#lineargrad = dwg.defs.add(svgwrite.gradients.LinearGradient(id="rwb"))
            g.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100),ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1, fill='url(#'+str(gradient_idc)+')'))
            #g.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100),ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1, fill='url(#'+str(gradient_idc)+')', mask='url(#mask1)'))
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


