from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys
from itertools import combinations
from random import sample
from svgwrite.container import Group
from Bio import SeqIO
from collections import defaultdict, deque
from scaffold import Scaffold, Longreads, LongReadSVG
import networkx as nx


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


scafs = Longreads(args.inputfiles, blacklist, args.linename, whitelist_lreads)
if whitelist_ctgs:
    scafs.filter_whitelist_ctgs(whitelist_ctgs)
scafs.filter_contigcounts(args.mincontigs)
reverse_mappers = set()
reverse_mappers.add("344DBB")
reverse_mappers.add("472DBB")
scafs.turn_longreads_around(reverse_mappers)
scafs.sort_by_starts()
scafs.filter_small_contigs(300)
#scafs.filter_overlapped_contigs(0.5)
scafs.filter_contigcounts(args.mincontigs) 


print("Reads meeting criteria: " + str(len(scafs.lreads)))
print("Pseudoaligning all...")
lr_scores, lr_dists = scafs.pseudoalign_all()
print("Pseudoaligning finished")

def N(node, dists):
    if node:
        return set(dists[node].keys()) - set([node])
    else:
        return set()

def BronKerbosch2(R,P,X,cliques,dists):
    if (not P) and (not X):
        cliques.append(R)
    if P:
        u = sample(P | X,1)[0]
    else:
        u = None
    for v in (P - N(u,dists)):
        Rn = R | set([v])
        Pn = (P-Rn) & N(v, dists)
        Xn = ((X & N(v, dists)) - Rn) - Pn
        BronKerbosch2(Rn, Pn, Xn,cliques, dists)
        P -= set(v)
        X |= set(v)

all_nodes = set(scafs.lreads.keys())

#print("Identifying cliques...")
#creads = []
#while all_nodes:
    #cliques = []
    #BronKerbosch2(set(),all_nodes,set(),cliques, lr_dists)
    #scliques = sorted(cliques, key = lambda x: len(x))
    #cluster = scliques[-1]
    #origin = sample(cluster,1)[0]
    #for read in cluster:
        #if lr_dists[read][origin] > 0:
            #origin = read
    #if len(cluster) == 1: 
        #break
    #creads.append(deque(sorted(cluster, key=lambda x: lr_dists[origin][x])))
    #all_nodes -= scliques[-1]
#print("Clique identification finished...")



# greedyly expand component with best matching thing
gr = nx.DiGraph()
for lr1i, lr1 in lr_scores.items():
    gr.add_edge(lr1i,lr1i,dist=0)
    if lr1i not in gr.nodes():
        gr.add_node(lr1i)
    ms = max(lr1.values())
    lr2i = max(lr1.keys(), key=lambda x: lr1[x])
    if lr2i not in gr.nodes():
        gr.add_node(lr2i)
    #print(lr2i)
    #print(lr_dists[lr1i][lr2i])
    if lr_dists[lr1i][lr2i] is not None:
        gr.add_edge(lr1i,lr2i,dist=lr_dists[lr1i][lr2i])
        gr.add_edge(lr2i,lr1i,dist=lr_dists[lr2i][lr1i])


# order components
creads = []
for component in list(nx.connected_components((gr.to_undirected()))):
    # get all distances from random anchor node
    if len(component) == 1:
        creads.append(deque(component))
        continue
    anchor = sample(component,1)[0]
    worklist = deque([anchor])
    ns = gr.neighbors(anchor)
    for n in ns:
        worklist.append(n)
    while worklist:
        #print("worklist: " + str(len(worklist)))
        cn = worklist.popleft()
        ns = nx.neighbors(gr,cn)
        #print("ns: " + str(len(ns)))
        for n in ns:
            if (anchor, n) not in gr.edges():
                gr.add_edge(anchor,n, dist=gr[anchor][cn]["dist"] + gr[cn][n]["dist"])
                gr.add_edge(n, anchor, dist= - gr[anchor][cn]["dist"] - gr[cn][n]["dist"])
                worklist.append(n)
    # now get all distances from the leftmost node
    minnode = min(gr[anchor], key = lambda x: gr[anchor][x]["dist"])
    for node in set(component):
        gr.add_edge(minnode, node, dist = gr[minnode][anchor]["dist"] + gr[anchor][node]["dist"])
        gr.add_edge(node, minnode, dist = -gr[minnode][anchor]["dist"] - gr[anchor][node]["dist"])
    # put this in data structure that's used by the plotting, update lr_dists
    creads.append(deque(sorted(component, key=lambda x: gr[minnode][x]["dist"])))
    #print(g[minnode])

            


print("Number of clusters found: " + str(len(creads)))
for item in creads:
    print("Number of reads: " + str(len(item)))

#sorted_reads = []
#smids = []
#smoffs = []
#sorted_clusters = []
#if args.alignreads:
    #for cluster in creads:
        #clst = cluster.copy()
        #origin = clst.popleft()
        #distance_known = set([origin])
        #while clst:
            #citem = clst.popleft()
            #if citem not in lr_dists[origin]:
                #for bitem in distance_known:
                    #if bitem in lr_dists[citem]:
                        #lr_dists[origin][citem] = lr_dists[origin][bitem] + lr_dists[bitem][citem]
                        #lr_dists[citem][origin] = -(lr_dists[origin][bitem] + lr_dists[bitem][citem])
            #distance_known.add(citem)
        ##print(lr_dists[origin])
        ##print(len(lr_dists[origin]))
        #sorted_reads = sorted(set(lr_dists[origin].keys()) & set(cluster), key = lambda x: lr_dists[origin][x])
        #smid = sorted_reads[0]
        ##print(smid)
        #for item in lr_dists[origin]:
            #if item not in lr_dists[smid]:
                #lr_dists[smid][item] = lr_dists[smid][origin] + lr_dists[origin][item]
        ##print(lr_dists[smid])
        ##print(len(lr_dists[smid]))
#
        #sorted_clusters.append(sorted_reads)



# draw interesting reads
ypos = 0
xtext = 10
xpad = 200
ypad = 10

image = LongReadSVG(args.output, zoom=100)
dwg = image.dwg

def shortname(ctgname):
    if "_" in ctgname:
        return ctgname.split("_")[1]
    else:
        return ctgname

csize = 0
sorted_clusters = creads

gradient_idc = 0

ypos += 200
for cluster in sorted_clusters:
    for rid in cluster:
        if args.alignreads:
            xoffset = gr[cluster[0]][rid]["dist"]
        else:
            xoffset = 0
        ypos += 28
        dwg.add(dwg.text(rid, insert=(xtext, ypad+ypos+1), fill='black', style="font-size:7"))
        g = dwg.defs.add(Group(id=rid))
        g.add(dwg.line((xpad+ xoffset/100, ypad+ypos), ( xpad + (xoffset+scafs.lreads[rid]["length"])/100, ypad+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))
        above = True
        col = "black"
        for ctg in scafs.lreads[rid]["maps"]:
            sc = ctg["scr"]
            ec = ctg["ecr"]
            scc = ctg["scc"]
            ecc = ctg["ecc"]
            ctgnr = ctg["name"].rstrip("rc").rstrip(args.linename)
            ctgn = ctg["name"].rstrip("rc")
            if ctgn.startswith("chr"):
                ctgn = ctgn[0:8]
                continue
            lenc = contigs[shortname(ctgn)]

            gradient_idc += 1
            #clip_path = dwg.defs.add(dwg.clipPath(id=str(clip_path_idc)))
            #clip_path.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100), ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1))
            #leftclip = scc/100
            #rightclip = (contigs[ctgn]-ecc)/100
            xa = -scc/(ecc-scc)
            xe = 1+(contigs[shortname(ctgn)]-ecc)/(ecc-scc)
            if ctg["strand"] == 1:
                xa = -(contigs[shortname(ctgn)]-ecc)/(ecc-scc)
                xe = 1+scc/(ecc-scc)
            lineargrad = dwg.defs.add(svgwrite.gradients.LinearGradient(id=str(gradient_idc), x1=xa, x2=xe, y1=0, y2=0))
            if ctgnr.endswith("0") or ctgnr.endswith("1"):
                col1 = "#FF0000"
                col2 = "#0000FF"
            elif ctgnr.endswith("2") or ctgnr.endswith("3"):
                col1 = "#FFE119"
                col2 = "#000000"
            elif ctgnr.endswith("4") or ctgnr.endswith("5"):
                col1 = "#911eb4"
                col2 = "#a9a9a9"
            elif ctgnr.endswith("6") or ctgnr.endswith("7"):
                col1 = "#000075"
                col2 = "#f58231"
            elif ctgnr.endswith("8") or ctgnr.endswith("9"):
                col1 = "#e6beff"
                col2 = "#808000"
                
            if ctg["strand"] == 1:
                tcol = col1
                col1 = col2
                col2 = tcol
            upperl = 0.90
            lowerl = 0.70
            opacity = ctg["quality"]/(upperl-lowerl) - lowerl/(upperl-lowerl)
            #opacity = 0.2
            opacity = min(max(opacity,0),1)
            lineargrad.add_stop_color("0%",col1,opacity)
            lineargrad.add_stop_color("50%","#FFFFFF",opacity)
            lineargrad.add_stop_color("100%",col2, opacity)
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
            g.add(dwg.text("s"+ctgn+"$", insert=(xpad+(xoffset+sc)/100,yt),fill=col, style="font-size:4"))
            if ctg["strand"] == 0:
                direction = ">"
            elif ctg["strand"] == 1:
                direction = "<"
            g.add(dwg.text(direction, insert=(xpad+(xoffset+sc)/100,ypad+ypos+2),style="font-size:6"))
        dwg.add(g)
    ypos += 28
    dwg.add(dwg.line((10, ypad+ypos), ( 10000, ypad+ypos+10), stroke=svgwrite.rgb(0, 0, 0, '%')))


dwg.save()


