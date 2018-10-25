from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys



parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("--whitelist", help="Only plot long reads with ids found in this whitelist file.")
parser.add_argument("--blacklist", help="Do not plot long reads in this file.")
parser.add_argument("--alignreads", action="store_true", help="Reads will be turned around if necessary and aligned according to their distances.")
parser.add_argument("cellline", help="Name of the cellline")
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
        if item["contig"].endswith(args.cellline):
            counter += 1
    if counter >= 2:
        greadst[rid] = read
        for item in read["maps"]:
            if item["contig"].endswith(args.cellline):
                greads_contigset.add(item["contig"])
    if counter == 1:
        intreads[rid] = read
        for item in read["maps"]:
            if item["contig"].endswith(args.cellline):
                intreads_contigset.add(item["contig"])
            
intersection = greads_contigset.intersection(intreads_contigset)
sintersection = sorted(intersection, key = lambda x: int(x.rstrip(args.cellline)))
#print(sintersection)


for rid in intreads:
    if len(intreads[rid]["maps"]) == 1:
        pass
        #print("\t".join([rid, str(intreads[rid])]))
        

# sort contigs by left coordinate
for rid in greadst:
    soverlaps = sorted(lreads[rid]["maps"], key = itemgetter("scr"))
    greads[rid] = greadst[rid]
    greads[rid]["maps"]=soverlaps

ogreads = greads.copy()
print("lgnth: " + str(len(greadst)))


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


# draw interesting reads
ypos = 0
xtext = 10
xpad = 200
ypad = 10
dwg = svgwrite.Drawing(args.output, size=(u'1700', u'4600'), profile='full')
csize = 0
bold = False
for cluster in creads:
    bold = not bold
    for rid in creads[cluster].items():
        ypos += 20
        if bold:
            dwg.add(dwg.text(rid[0], insert=(xtext, ypad+ypos+1), fill='red', style="font-size:7"))
        else:
            dwg.add(dwg.text(rid[0], insert=(xtext, ypad+ypos+1), fill='black', style="font-size:7"))
        dwg.add(dwg.line((xpad, ypad+ypos), ( xpad + ogreads[rid[0]]["length"]/100, ypad+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))
        above = True
        col = "black"
        for read in ogreads[rid[0]]["maps"]:
            #print(read)
            sc = read["scr"]
            ec = read["ecr"]
            ctg = read["contig"].rstrip(args.cellline)
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctg = ctg[0:8]
            dwg.add(svgwrite.shapes.Rect((xpad+(sc/100),ypad+ypos-3), ((ec-sc)/100,6), stroke='black', stroke_width=1, fill = 'white'))
            if above:
                yt = ypad+ypos-4
                if col == "blue":
                    col = "black"
                else: 
                    col = "blue"
            else:
                yt = ypad+ypos+7
            above = not above
            dwg.add(dwg.text(ctg, insert=(xpad+sc/100,yt),fill=col, style="font-size:4"))
            if read["strand"] == 0:
                direction = ">"
            elif read["strand"] == 1:
                direction = "<"
            dwg.add(dwg.text(direction, insert=(xpad+sc/100,ypad+ypos+2),style="font-size:6"))
    dwg.add(dwg.line((10, ypad+ypos+10), ( 1000, ypad+ypos+10), stroke=svgwrite.rgb(0, 0, 0, '%')))


dwg.save()





