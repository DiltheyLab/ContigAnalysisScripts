import svgwrite
from svgwrite.container import Group
import sys
from Bio import SeqIO
from scaffold import LongReadSVG
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputfile", help="FASTA input file")
parser.add_argument("output", help="SVG output file")
args = parser.parse_args()

seqs = {}
for read in SeqIO.parse(args.inputfile, "fasta"):
    seqs[read.id] = str(read.seq).upper()

image = LongReadSVG(args.output, zoom=800)
dwg = image.dwg

ypad = 7
xpad = 20
col1 = "black"
col2 = "grey"
col = col1
coln = "red"
#ypos = yoff
y_size = 100
ypos = 200
ypos += y_size + ypad
y_halfsize = y_size/2
xoffset = 0 #no alignment

for id, seq in seqs.items():
    dwg.add(dwg.text(id, insert=(xpad - 500, ypad+ypos+1), fill='black', style="font-size:40"))
    #dwg.add(dwg.text(scafs.lreads[rid]["fname"], insert=(xtext, ypad+ypos+5), fill='black', style="font-size:4"))
    g = dwg.defs.add(Group(id=id))
    status = "n" if seq[0] == "N" else "s"
    sidx = 0
    idx = 0
    print("id: " + str(id))
    while idx < len(seq):
        if status == "s":
            while idx < len(seq):
                if seq[idx] != "N":
                    idx += 1
                    continue
                else:
                    g.add(svgwrite.shapes.Rect((xpad+((xoffset+sidx)/image.zoom),ypad+ypos-y_halfsize), ((idx-1-sidx)/image.zoom,y_size), stroke='black', stroke_width=1, fill=col))
                    lcol = col
                    col = "red"
                    sidx = idx
                    status = "n"
                    idx += 1
                    break
                if idx == len(seq):
                    break
        else:
            while idx < len(seq):
                if seq[idx] == "N":
                    idx += 1
                    continue
                else:
                    g.add(svgwrite.shapes.Rect((xpad+((xoffset+sidx)/image.zoom),ypad+ypos-y_halfsize), ((idx-1-sidx)/image.zoom,y_size), stroke='black', stroke_width=1, fill=col))
                    col = col1 if lcol == col2 else col2
                    sidx = idx
                    status = "s"
                    idx += 1
                    break
                idx += 1
                if idx == len(seq):
                    break

    dwg.add(g)
    ypos += 2*y_size

dwg.save()
