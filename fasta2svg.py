import svgwrite
from svgwrite.container import Group
import sys
import re
from Bio import SeqIO
from scaffold import LongReadSVG, Longreads
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputfile", help="FASTA/PAF input file")
parser.add_argument("output", help="SVG output file")
parser.add_argument("--switch_colors", action="store_true", default=False)
args = parser.parse_args()

altnames = {"chr6_GL000250v2_alt": "APD",
            "chr6_GL000251v2_alt": "COX",
            "chr6_GL000252v2_alt": "DBB",
            "chr6_GL000253v2_alt": "MANN",
            "chr6_GL000254v2_alt": "MCF",
            "chr6_GL000255v2_alt": "QBL",
            "chr6_GL000256v2_alt": "SSTO"}

seqs = {}
format = "fasta"

with open(args.inputfile) as f:
    for line in f:
        if line.startswith(">"):
            break
        else:
            format = "paf"
            celllinestr = line.split()[0]
            m = re.search('[A-Z]+', celllinestr)
            cellline = m.group(0)
            print("cell line detected: " + cellline)
            break

if format == "fasta":
    for read in SeqIO.parse(args.inputfile, "fasta"):
        seqs[read.id] = str(read.seq).upper()
elif format == "paf":
    scaf = Longreads.init_from_reverse_paf(args.inputfile)
    scaf.sort_by_starts()
    lread = scaf.lreads

else:
    print("Problem! Format unknown")
    sys.exit()
    

image = LongReadSVG(args.output, zoom=800)
dwg = image.dwg

ypad = 7
xpad = 20
col1 = "black"
col2 = "lightgrey"
col = col1
coln = "red"
#ypos = yoff
y_size = 200
ypos = 70
ypos += y_size + ypad
y_halfsize = y_size/2
xoffset = 0 #no alignment
cat_contigs = True
last_ecr = 0

if format == "fasta":
    for id, seq in seqs.items():
        totalNs = 0
        if id in altnames:
            dwg.add(dwg.text(altnames[id], insert=(xpad - 500, ypad+ypos+1), fill='black', style="font-size:40"))
        else:
            dwg.add(dwg.text(id, insert=(xpad - 500, ypad+ypos+1), fill='black', style="font-size:40"))
        #dwg.add(dwg.text(scafs.lreads[rid]["fname"], insert=(xtext, ypad+ypos+5), fill='black', style="font-size:4"))
        g = dwg.defs.add(Group(id=id))
        status = "n" if seq[0] == "N" else "s"
        sidx = 0
        idx = 0
        #print("id: " + str(id) + "\t" + str(altnames[id]))
        while idx < len(seq):
            if status == "s":
                while idx < len(seq):
                    if seq[idx] != "N":
                        idx += 1
                        continue
                    else:
                        g.add(svgwrite.shapes.Rect((xpad+((xoffset+sidx)/image.zoom),ypad+ypos-y_halfsize), ((idx-1-sidx)/image.zoom,y_size), stroke='black', stroke_width=0, fill=col))
                        sidx = idx
                        status = "n"
                        idx += 1
                        break
                    if idx == len(seq):
                        break
            elif status == "n":
                while idx < len(seq):
                    if seq[idx] == "N":
                        idx += 1
                        continue
                    else:
                        g.add(svgwrite.shapes.Rect((xpad+((xoffset+sidx)/image.zoom),ypad+ypos-y_halfsize), ((idx-1-sidx)/image.zoom,y_size), stroke='black', stroke_width=0, fill="red"))
                        col = col1 if col == col2 else col2
                        totalNs += (idx - sidx)
                        sidx = idx
                        status = "s"
                        idx += 1
                        break
                    idx += 1
                    if idx == len(seq):
                        break
        if seq[idx-1] == "N":
            totalNs += (idx -sidx)
        dwg.add(g)
        print("number of Ns: " + str(totalNs))
        ypos += 2*y_size
elif format == "paf":
    for rid, lr in lread.items():
        print(rid)
        dwg.add(dwg.text(rid, insert=(xpad - 300, ypad+ypos+1), fill='black', style="font-size:80"))
        g = dwg.defs.add(Group(id=rid))
        #red rectangle background
        lc = lr["maps"][-1]["ecr"]
        g.add(svgwrite.shapes.Rect((xpad+((xoffset)/image.zoom),ypad+ypos-y_halfsize), ((lc)/image.zoom,y_size), stroke='black', stroke_width=0, fill=coln))
        totalNs = 0
        for ctg in lr["maps"]:
            scr = ctg["scr"]
            length = ctg["ecr"] -scr
            g.add(svgwrite.shapes.Rect((xpad+((xoffset+scr)/image.zoom),ypad+ypos-y_halfsize), ((length)/image.zoom,y_size), stroke='black', stroke_width=0, fill=col))
            #totalNs -= length
            if ctg != lr["maps"][0] and last_ecr > ctg["scr"]:
                totalNs += last_ecr - ctg["scr"]

            if cat_contigs:
                if last_ecr < ctg["scr"]:
                    col = col1 if col == col2 else col2
            else:
                col = col1 if col == col2 else col2
            last_ecr = ctg["ecr"] 
        dwg.add(g)
        print("number of Ns: " + str(totalNs))
        ypos += 2*y_size
        

dwg.save()
