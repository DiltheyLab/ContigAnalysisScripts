from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import random
from operator import itemgetter
from itertools import combinations, cycle, product
from collections import defaultdict
import svgwrite
import squarify
from scaffold import Scaffold, Scaffolds

parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("--paf", help="Input is paf file", action="store_true", default = False)
parser.add_argument("--mincontigs", type=int, default=2,help="Minimum number of contigs on long read for the read to be considered")
parser.add_argument("--summaryfile", help="Contig Distance Summary file")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
parser.add_argument("--mergefile", help="File that contains merging information to retrieve the sequence of scaffolds after.")
parser.add_argument("--contigsmergefile", help="File with information on what contigs should be merged.")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("linename", help="Name of cell line")
parser.add_argument("SVG", help="Scaffolds are drawn to this SVG file")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=25)

args = parser.parse_args()

reads = {}
greads = {}
contigs = {}
allcontigs = {}
contig2scaffold = defaultdict(list)
cluster_counter = 1

srneighs = dict()

if args.mergefile:
    with open(args.mergefile, "w+") as mergef:
       pass 

ctgpos = {}
pos = 0
for read in SeqIO.parse(args.contigfile, "fasta"):
    pos += 1
    contigs[read.id] = len(read.seq)
    allcontigs[read.id] = len(read.seq)
    ctgpos[read.id] = pos


blacklist = defaultdict(list)
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            if sline[0] == "contig":
                blacklist[sline[1]] = "y"
            else:
                blacklist[sline[0]].append(sline[1])
#print(blacklist)

def get_other_relpos(relpos):
    if relpos == "left":
        return "right"
    return "left"

def add_neighs(ctg1, ctg2, relpos, dist):
    if ctg1 in srneighs:
        if ctg2 not in [x[0] for x in srneighs[ctg1][relpos]]:
            srneighs[ctg1][relpos].append((ctg2,dist))
        if ctg2 in [x[0] for x in srneighs[ctg1][get_other_relpos(relpos)]]:
            pass
            #print("Conflict for " + ctg1 + " and " + ctg2)
            #print(srneighs[ctg1])
    else:
        srneighs[ctg1] = {"left": [], "right": []}
        srneighs[ctg1][relpos] = [(ctg2, dist)]

if args.summaryfile:
    with open(args.summaryfile) as f:
        for line in f:
            sline = line.split()
            #if sline[1] == "NA":
            #    continue
            if int(sline[2]) < args.mindepth:
                continue
            ctg1 = sline[0].split("_")[0].strip("+").strip("-")
            ctg2 = sline[0].split("_")[1].strip("+").strip("-")
            if ctg1 in blacklist or ctg2 in blacklist:
                continue
            ori1 = sline[0].split("_")[0][0]
            ori2 = sline[0].split("_")[1][0]
            distance = float(sline[3])
            dist = int(distance)
            if distance < 0 and abs(distance) > allcontigs[ctg1]:
                continue
            if distance < 1000 and distance > -1000:
                if ori1 == ori2:
                    if ori1 == "+":
                        add_neighs(ctg1,ctg2,"right",dist)
                        add_neighs(ctg2,ctg1,"left",dist)
                    else:
                        add_neighs(ctg2,ctg1,"right",dist)
                        add_neighs(ctg1,ctg2,"left",dist)
                else: # some contigs are defined in wrong orientation
                      # we still want to handle these
                    if ori1 =="+": # so ori2 == "-"  |ctg1 >| |< ctg2|
                        add_neighs(ctg1,ctg2,"right",dist)
                        add_neighs(ctg2,ctg1,"right",dist)
                    else: # |ctg1 <| |> ctg2|
                        add_neighs(ctg1,ctg2,"left",dist)
                        add_neighs(ctg2,ctg1,"left",dist)

# get rid of prefixes a_ b_ c_ and the like
def shortname(ctgname):
    if "_" in ctgname:
        return ctgname.split("_")[1]
    else:
        return ctgname

# some custom additions
#add_neighs("224APD","976APD","right",round(16.67))
#add_neighs("916APD","1247APD","right",round(116.8))
#add_neighs("2406APD","1671APD","right",round(52.3))
#add_neighs("105APD","726APD","right",round(-40.1))
#add_neighs("726APD","1674APD","right",round(304.0))
#add_neighs("2080APD","2377APD","right",round(361.9))
#add_neighs("2377APD","928APD","right",round(75.6))
#add_neighs("504APD", "49APD","right",round(-460.3))
#add_neighs("49APD","116APD","right",round(-4.8))
#srneighs["504APD"]["right"] = []

   


print("Nr. of scaffolds: " + str(len(contigs)))

lrs = Scaffolds(args.efile, args.paf, blacklist, args.linename)
lrs.filter_contigcounts(args.mincontigs)
scaffolds = lrs.construct_scaffolds(allcontigs)

#scaffolds = {}
#for rid in scafs.lreads:
    #nscaff = Scaffold.init_from_LR((rid,scafs.lreads[rid]),args.linename, args.paf, allcontigs )
for scafid, scaf in scaffolds.items():
    for ctg in scaf.contigset:
        if ctg.endswith(args.linename):
            contig2scaffold[ctg].append(id(scaf))
        
    #nscaff.get_sequence_th()
    #nscaff.add_seq_info()
    #scaffolds[id(nscaff)] = nscaff


toRemove = set()
for idx,scaf in scaffolds.items():
    #find and get rid of conflicting scaffolds
    #if scaf.find_conflicts():
    #    toRemove.add(idx)
    octgs = scaf.remove_overlapped_contigs(allcontigs)
    sctgs = scaf.remove_short_contigs(allcontigs)
    for ctg in octgs | sctgs:
        contig2scaffold[ctg].remove(idx)
    if len(scaf.contigset) == 0:
        toRemove.add(idx)
        
    
for idx in toRemove:
    del(scaffolds[idx])
    



dwg = svgwrite.Drawing(args.SVG,size=(u'1700', u'2200'), profile='full')

lineargrad = dwg.defs.add(svgwrite.gradients.LinearGradient(id="rwb"))
lineargrad.add_stop_color(0,"#FF0000")
lineargrad.add_stop_color(50,"#FFFFFF")
lineargrad.add_stop_color(100,"#0000FF")
#svgwrite.gradients.LinearGradient((0,0), (100,0))
yp = 10
xtext = 10
xpad = 20
dwg.add(dwg.text("600 bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 5 
dwg.add(dwg.line((xpad, yp), ( xpad + 600/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad + 600/100, yp-2), ( xpad +600/100, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
yp += 10
dwg.add(dwg.text("10000 bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 5
dwg.add(dwg.line((xpad, yp), ( xpad + 10000/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
for i in range(1,10):
    dwg.add(dwg.line( (xpad + (10000/100)/10 * i, yp-1), (xpad + (10000/100)/10 * i, yp+1), stroke=svgwrite.rgb(0,0,0,'%')))
dwg.add(dwg.line((xpad + 10000/100, yp-2), ( xpad +10000/100, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
yp += 10
dwg.add(dwg.text("100k bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 5
dwg.add(dwg.line((xpad, yp), ( xpad + 100000/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
for i in range(1,10):
    dwg.add(dwg.line( (xpad + (100000/100)/10 * i, yp-1), (xpad + (100000/100)/10 * i, yp+1), stroke=svgwrite.rgb(0,0,0,'%')))
dwg.add(dwg.line((xpad + 100000/100, yp-2), ( xpad +100000/100, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
yp += 10
g1M = dwg.defs.add(dwg.g(id='g003'))
g1M.add(dwg.text("M bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 5
g1M.add(dwg.line((xpad, yp), ( xpad + 3000000/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%'), stroke_width="5"))
g1M.add(dwg.line((xpad, yp-7), ( xpad , yp+7), stroke=svgwrite.rgb(0, 0, 0, '%'), stroke_width="5"))
for i in range(0,31):
    if i%10 == 0:
        g1M.add(dwg.line( (xpad + (1000000/100)/10 * i, yp-7), (xpad + (1000000/100)/10 * i, yp+200), stroke=svgwrite.rgb(0,0,0,'%'),stroke_width="5"))
    else:
        g1M.add(dwg.line( (xpad + (1000000/100)/10 * i, yp-5), (xpad + (1000000/100)/10 * i, yp+200), stroke=svgwrite.rgb(0,0,0,'%'),stroke_width="2"))
g1M.add(dwg.text("1 M", insert=( xpad+ 1000000/100-10, yp-10), fill='black', style="font-size:10"))
g1M.add(dwg.text("2 M", insert=( xpad+ 2000000/100-10, yp-10), fill='black', style="font-size:10"))
g1M.add(dwg.text("3 M", insert=( xpad+ 3000000/100-10, yp-10), fill='black', style="font-size:10"))
dwg.add(g1M)
yp += 20
rect = dwg.add(svgwrite.shapes.Rect((xpad,yp-3), (2000/100,7), stroke='green', stroke_width=1 ))
rect.fill(color="none").dasharray([2, 2])
dwg.add(dwg.text("region of scaffold covered by long read", insert=( xpad+ 2000/100 + 5, yp+2), fill='black', style="font-size:7"))
yp += 10
if args.summaryfile:
    dwg.add(svgwrite.shapes.Rect((xpad,yp-3), (2000/100,7), stroke='gray', stroke_width=1, fill='none' ))
    dwg.add(dwg.text("contigs from short read data", insert=( xpad+ 2000/100 + 5, yp+2), fill='black', style="font-size:7"))
    yp += 10
dwg.add(svgwrite.shapes.Rect((xpad,yp-3), (2000/100,7), stroke='black', stroke_width=1, fill='none' ))
dwg.add(dwg.text("contigs from long read data", insert=( xpad+ 2000/100 + 5, yp+2), fill='black', style="font-size:7"))
yp += 20

# cluster np-reads 
print("scaffolding long reads ....")
creads = {}
clusternr = 0
olen_scaf = len(scaffolds)+1
while len(scaffolds)-olen_scaf != 0:
    olen_scaf = len(scaffolds)
    for contig in contig2scaffold:
        if len(contig2scaffold[contig]) > 1:
            scaf1 = scaffolds[contig2scaffold[contig][0]]
            if contig2scaffold[contig][1] not in scaffolds:
                print("here is the problem: " + str(ctg) + " is not in c2s")
            scaf2 = scaffolds[contig2scaffold[contig][1]]
            success = scaf1.merge(scaf2, args.mergefile)
            for ctg in scaf2.contigset:
                if success:
                    if id(scaf1) not in contig2scaffold[ctg]:
                        contig2scaffold[ctg].append(id(scaf1))
                if id(scaf2) in contig2scaffold[ctg]:
                    contig2scaffold[ctg].remove(id(scaf2))
            del(scaffolds[id(scaf2)])
            break

print("Nr. of scaffolds: " + str(len(scaffolds)))
# add merge info to mergefile 
#if args.mergefile:
#    for scafid, scaf in scaffolds.items():
#        if not scaf.in_mergefile:
#            mode = "not_merged"
#            with open(args.mergefile, "a+") as mergef:
#                mergef.write("\t".join([mode ,scaf.name, str(scaf.turned_around)]))
#                mergef.write("\n")
#        print("Scaffold " + str(scaf.name))
#        print(", ".join(list(scaf.contigset)))


def stringify(nr):
    if nr >= 1000000:
        return "{:.1f}".format(nr/1000000) + " M"
    elif nr >= 10000:
        return "{:.0f}".format(nr/1000) + " k"
    elif nr > 1000:
        return "{:.1f}".format(nr/1000) + " k"
    else:
        return str(nr)

lr_lengths = []
lr_lengths.append([])
lr_lengths.append([])
for scaf in scaffolds.values():
    lr_lengths[1].append(scaf.length)
    #lr_lengths[0].append(scaf.name + "\n" + stringify(scaf.length))
    if scaf.length > 10000:
        lr_lengths[0].append(stringify(scaf.length))
    else:
        lr_lengths[0].append("")
    #lr_lengths[0].append(scaf.name)
    #print("length of scaffold " + scaf.name + ": " + str(scaf.length))
norm = matplotlib.colors.Normalize(vmin=min(lr_lengths[1]), vmax=max(lr_lengths[1]))
colors = [matplotlib.cm.gist_ncar(value) for value in random.sample(range(1024),len(lr_lengths[1]))]
lr_lengths.append(colors)


def is_rightmost(ctg):
    try:
        assert(len(contig2scaffold[ctg]) == 1)
    except AssertionError:
        print("Other than one scaffold for contig " + ctg1 + " | number of contigs: " + str(len(contig2scaffold[ctg])))
    scaf = scaffolds[contig2scaffold[ctg][0]]
    if ctg == scaf.get_rightmost_contig():
        return True
    return False

def is_leftmost(ctg):
    try:
        assert(len(contig2scaffold[ctg]) == 1)
    except AssertionError:
        print("Other than one scaffold for contig " + ctg1 + " | number of contigs: " + str(len(contig2scaffold[ctg])))
    #try:
    #print(contig2scaffold[ctg])
    scaf = scaffolds[contig2scaffold[ctg][0]]
    #except KeyError:
    #    print("Problem during is_leftmost. For contig: " + ctg)
    #    sys.exit(1)
    if ctg == scaf.get_leftmost_contig():
        return True
    return False

#sys.exit()

if args.summaryfile:
    print("adding short reads ....")
    change = True
    while change:
        change = False
        for scafid in scaffolds.keys():
            scaffold = scaffolds[scafid]
            ctg1 = scaffold.get_rightmost_contig()
            if not ctg1 in srneighs:
                continue
            for ctg2n,ctg2d in srneighs[ctg1]["right"]:
                if ctg2n in contig2scaffold and len(contig2scaffold[ctg2n]) > 0 and is_leftmost(ctg2n):
                #if ctg2n in contig2scaffold and len(contig2scaffold[ctg2n]) > 0 :
                    scaf2 = scaffolds[contig2scaffold[ctg2n][0]]
                    #print("merged away " + str(id(scaf2)))
                    #print("contig: " + ctg1)
                    scaf1 = scaffolds[contig2scaffold[ctg1][0]]
                    if scaf1.name == scaf2.name:
                        continue
                    if scaffold.merge_sr(ctg1,ctg2n,scaf1, scaf2, ctg2d,allcontigs, args.mergefile):
                        del(scaffolds[id(scaf2)])
                        for ctg in scaf2.contigset:
                            contig2scaffold[ctg] = [id(scaf1)]
                            scaf1.contigset.add(ctg)
                        for ctg in scaf2.contigset_sr:
                            contig2scaffold[ctg] = [id(scaf1)]
                            scaf1.contigset_sr.add(ctg)
        
                    
                    change = True
                    break
                elif ctg2n in contig2scaffold and len(contig2scaffold[ctg2n]) > 0 :
                    #print("Probably " + ctg1 + " fits on " + scaf2.name + " left of " + ctg2n)
                    pass
                else:
                    scaffold.add_short_read_contig_right(ctg1,ctg2n,ctg2d, 0, allcontigs, args.mergefile)
                    contig2scaffold[ctg2n].append(id(scaffold))
                    if ctg2n not in contigs:
                        print("Contig " + ctg2n + " is already part of a scaffold. Investigate!")
                    else:
                        del(contigs[ctg2n])
                    change = True
            #for ctg2n,ctg2d in srneighs[ctg1]["left"]:
            #    if ctg2n in contig2scaffold:
            #        pass
            #        #print("Problem merging")
            #    else:
            #        scaffold.add_short_read_contig_left(ctg1,ctg2n,ctg2d, 0, args.mergefile)
            #        if ctg2n not in contigs:
            #            print("Contig " + ctg2n + " is already part of a scaffold. Investigate!")
            #        else:
            #            del(contigs[ctg2n])
            #        change = True
            if change:
                break

print("Nr. of scaffolds: " + str(len(scaffolds) + len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")
        
if args.contigsmergefile:
    with open(args.contigsmergefile, "w+") as cmergef:
        for scid, scaffold in scaffolds.items():
            sortedcontigs = sorted(scaffold.contigset | scaffold.contigset_sr, key = lambda item: scaffold.left_coords[item])
            cmergef.write(">" + scaffold.name + "\n")
            for ctg in scaffold.contigset:
                if ctg.startswith("chr"):
                    sortedcontigs.remove(ctg)
            for ctg1, ctg2 in zip(sortedcontigs[:-1], sortedcontigs[1:]):
                cmergef.write(ctg1+ "\t" + ctg2 + "\n")
            
# Draw all scaffolds
lrsr_lengths = []
lrsr_lengths.append([])
lrsr_lengths.append([])
for scaf in scaffolds.values():
    lrsr_lengths[1].append(scaf.length)
    lrsr_lengths[0].append(stringify(scaf.length))
    #lrsr_lengths[0].append("c_" + scaf.name.split("_")[1])
    yp += scaf.to_SVG(dwg, allcontigs,  xpad, yp, False) + 10
dwg.save()
#norm = matplotlib.colors.Normalize(vmin=min(lrsr_lengths[1]), vmax=max(lrsr_lengths[1]))
#colors = [matplotlib.cm.tab20b(norm(value)) for value in lrsr_lengths[1]]
#lrsr_lengths.append(colors)
plt.rc('font', size=15)          # controls default text sizes
plt.subplot(121)
squarify.plot(sizes=lr_lengths[1], label=lr_lengths[0], alpha=.9,color = colors )
#squarify.plot(sizes=lr_lengths[1], alpha=.9 )
plt.axis('off')
plt.title('after long read scaffolding')
plt.subplot(122)
squarify.plot(sizes=lrsr_lengths[1], label=lrsr_lengths[0], alpha=.9 )
#squarify.plot(sizes=lrsr_lengths[1], alpha=.9 )
plt.axis('off')
plt.title('after long + short read scaffolding')
plt.show()

