from argparse import ArgumentParser
import svgwrite
from operator import itemgetter
import sys


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("--whitelist", help="Only plot long reads with ids found in this whitelist file.")
parser.add_argument("cellline", help="Name of the cellline")
parser.add_argument("output", help="SVG output file")

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
        strand = sline[8]
        sc = int(sline[5])
        ec = int(sline[6])
        sc_c = int(sline[9])
        ec_c = int(sline[10])
        c_l = int(sline[11])
        if rid in reads:
            reads[rid]["overlaps"].append([ctg,strand,sc,ec,sc_c,ec_c,c_l])
        else:
            reads[rid] = {}
            reads[rid]["length"] = int(sline[7])
            reads[rid]["overlaps"] = [[ctg,strand,sc,ec,sc_c,ec_c,c_l]]

# get interesting reads
greadst = {}
intreads = {}
greads_contigset = set()
intreads_contigset = set()

for rid in reads:
    counter = 0
    for item in reads[rid]["overlaps"]:
        if item[0].endswith(args.cellline):
            counter += 1
    if counter >= 2:
        greadst[rid] = reads[rid]
        for item in reads[rid]["overlaps"]:
            if item[0].endswith(args.cellline):
                greads_contigset.add(item[0])
    if counter == 1:
        intreads[rid] = reads[rid]
        for item in reads[rid]["overlaps"]:
            if item[0].endswith(args.cellline):
                intreads_contigset.add(item[0])
            
intersection = greads_contigset.intersection(intreads_contigset)
sintersection = sorted(intersection, key = lambda x: int(x.rstrip(args.cellline)))
#print(sintersection)


for rid in intreads:
    if len(intreads[rid]["overlaps"]) == 1:
        pass
        #print("\t".join([rid, str(intreads[rid])]))
        

# sort contigs by left coordinate
for rid in greadst:
    soverlaps = sorted(reads[rid]["overlaps"], key = itemgetter(2))
    #print(reads[rid]["overlaps"])
    greads[rid] = greadst[rid]
    greads[rid]["overlaps"]=soverlaps

ogreads = greads.copy()


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
        for mapping in read["overlaps"]:
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
            for contig in cr[1]["overlaps"]:
                if not contig[0].startswith("chr"):
                    current_contigs.add(contig[0])
            #print("contigs: " + str(current_contigs))
            for readid,readval in greads.items():
                contig_found = False
                for contig in readval["overlaps"]:
                    if contig[0] in current_contigs:
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
        for read in ogreads[rid[0]]["overlaps"]:
            #print(read)
            sc = read[2]
            ec = read[3]
            ctg = read[0].rstrip(args.cellline)
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
            if read[1] == "0":
                direction = ">"
            elif read[1] == "1":
                direction = "<"
            dwg.add(dwg.text(direction, insert=(xpad+sc/100,ypad+ypos+2),style="font-size:6"))
    dwg.add(dwg.line((10, ypad+ypos+10), ( 1000, ypad+ypos+10), stroke=svgwrite.rgb(0, 0, 0, '%')))


dwg.save()





