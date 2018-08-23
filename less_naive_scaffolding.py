from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from itertools import combinations
import svgwrite


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("--mincontigs", type=int, default=2,help="Minimum number of contigs on long read for the read to be considered")
parser.add_argument("--summaryfile", help="Contig Distance Summary file")
parser.add_argument("--blacklistfile", help="File containing long read ids where certain contig mappings should be ignored.")
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
contig2scaffold = {}


srneighs = dict()

def get_other_relpos(relpos):
    if relpos == "left":
        return "right"
    return "left"

def add_neighs(ctg1, ctg2, relpos, dist):
    if ctg1 in srneighs:
        if ctg2 not in [x[0] for x in srneighs[ctg1][relpos]]:
            srneighs[ctg1][relpos].append((ctg2,dist))
        if ctg2 in [x[0] for x in srneighs[ctg1][get_other_relpos(relpos)]]:
            print("Conflict for " + ctg1 + " and " + ctg2)
            print(srneighs[ctg1])
    else:
        srneighs[ctg1] = {"left": [], "right": []}
        srneighs[ctg1][relpos] = [(ctg2, dist)]


if args.summaryfile:
    with open(args.summaryfile) as f:
        for line in f:
            sline = line.split()
            if sline[1] == "NA":
                continue
            if int(sline[2]) < args.mindepth:
                continue
            ctg1 = sline[0].split("_")[0].strip("+").strip("-")
            ctg2 = sline[0].split("_")[1].strip("+").strip("-")
            ori1 = sline[0].split("_")[0][0]
            ori2 = sline[0].split("_")[1][0]
            distance = float(sline[1])
            dist = int(distance)
            if distance < 1000 and distance > -200:
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

blacklist = {}
if args.blacklistfile:
    with open(args.blacklistfile) as f:
        for line in f:
            sline = line.split()
            blacklist[sline[0]] = sline[1]

print(blacklist)
   
ctgpos = {}
pos = 0
for read in SeqIO.parse(args.contigfile, "fasta"):
    pos += 1
    contigs[read.id] = len(read.seq)
    allcontigs[read.id] = len(read.seq)
    ctgpos[read.id] = pos
    

print("Nr. of scaffolds: " + str(len(contigs)))

class Scaffold:
    scaf_info = []
    sr_info = dict()
    lr_info = dict()
    longread_coords = dict()
    coords = []
    left_coords = dict()
    right_coords = dict()
    orientation = dict()
    contigset = set()
    contigset_sr = set()
    length = 0
    sequence = ""
    nr_of_scaffolds = 0
    turned_around = False
    idx = ""
    get_sequence_fragments = []
    #id = 
    # All contigs will have scaffold coordinates.
    # Before scaffolds are merged all scaffold coordinates 
    # are equal to long read coordinates of the orginal long read 
    # that initialized the scaffold
    
    def __init__(self):
        Scaffold.nr_of_scaffolds += 1
        self.contigset = set()
        self.contigset_sr = set()
        #self.turned_around = False
        #lr_info[lr[0]] = lr[1]

    
    def delete(self):
        for ctg in self.contigset:
            contig2scaffold[ctg].remove(id(self))
            if contig2scaffold[ctg] == []:
                del(contig2scaffold[ctg])
            if not ctg in contigs:
                contigs[ctg] = allcontigs[ctg]
        del(scaffolds[self.idx])

    def print_contigset(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: int(item.rstrip(args.linename)))
        for contig in sortedcontigs:
            print(contig)

    def print_lrids(self):
        print(self.lr_info.keys())

    def print_coords(self):
        ostr = []
        for coord in self.coords:
            if coord != ("","",""):
                ostr.append(str(coord))
        print(" ".join(ostr))
        #print(self.coords)
        
    def print_contig_sequence(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: self.left_coords[item])
        print("-".join(sortedcontigs))
        
    def print_id(self):
        print(self.idx)

    # this signifies a hard turnaround
    # i.e. the coordinate system of the scaffold is reversed
    def turn_around(self):
        self.turned_around = True
        new_left_coords = dict()
        new_right_coords = dict()
        new_orientation = dict()
        
        for contig in self.contigset:
            new_right_coords[contig] = self.length - self.left_coords[contig] 
            new_left_coords[contig] = self.length - self.right_coords[contig] 
            new_orientation[contig] = 0 if self.orientation[contig] == 1 else 1
        self.left_coords = new_left_coords
        self.right_coords = new_right_coords
        self.orientation = new_orientation
        self.coords.reverse()

    def is_left_of(self,ctg1, ctg2):
        if self.left_coords(ctg1) < self.left_coords(ctg2):
            try:
                assert(self.right_coords(ctg1) < self.right_coords(ctg2))
            except AssertionError:
                print(ctg1 + " includes " + ctg2)
            return True
        else:
            try:
                assert(self.right_coords(ctg1) > self.right_coords(ctg2))
            except AssertionError:
                print(ctg2 + " includes " + ctg1)
            return False

    def get_leftmost_contig(self):
        ctgset = self.contigset.union(self.contigset_sr)
        sortedcontigs = sorted(ctgset, key = lambda item: self.left_coords[item])
        return(sortedcontigs[0])

    def get_rightmost_contig(self):
        ctgset = self.contigset.union(self.contigset_sr)
        sortedcontigs = sorted(ctgset, key = lambda item: self.right_coords[item])
        return(sortedcontigs[-1])


    def find_conflicts(self,max_dist):
        sortedcontigs = sorted(self.contigset, key = lambda item: self.left_coords[item])
        sortedcontigs2 = sorted(self.contigset, key = lambda item: self.right_coords[item])
        scstring1 = "-".join(sortedcontigs)
        scstring2 = "-".join(sortedcontigs2)
        if scstring1 != scstring2:
            print("Conflict! " + scstring1 + " " + scstring2 + " . Removing scaffold.")
            return True
        for ctg1,ctg2 in zip(sortedcontigs[:-1],sortedcontigs[1:]):
            if ctgpos[ctg1] > ctgpos[ctg2]:
                print("Contigs " + ctg1 + " and " + ctg2 + " ordered badly. Removing scaffold.")
                #return True
            dist = ctgpos[ctg2] - ctgpos[ctg1]
            if dist > max_dist:
                print("Contigs " + ctg1 + " and " + ctg2 + " too far apart (" + str(dist) + "). Removing scaffold.")
                #return True
        return False

    # Returns the space in y that it needs (depends on the number of longreads that were merged into this scaffold)
    def to_SVG(self, img, xoff, yoff):
        ypad = 7
        col = "black"
        nr_longreads = len(self.longread_coords)
        y_space_per_longread = 4
        ypos = 0
        for lrid, lrc in self.longread_coords.items():
            ylen = nr_longreads * y_space_per_longread - ypos
            rect = img.add(svgwrite.shapes.Rect((xoff+(lrc[0]/100),yoff+ypos), ((lrc[1]-lrc[0])/100,ylen+ypad), stroke='green', stroke_width=1 ))
            rect.fill(color="none").dasharray([2, 2])
            img.add(dwg.text(lrid, insert=(xoff+(lrc[0]/100),yoff+ypos-1),fill="green", style="font-size:2"))
            ypos += y_space_per_longread
        ypos += ypad
            
        img.add(dwg.line((xoff, yoff+ypos), ( xoff + self.length/100, yoff+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))

        above = True
        for ctg in sorted(self.contigset, key= lambda x: self.left_coords[x]):
            #print(read)
            sc = self.left_coords[ctg]
            ec = self.right_coords[ctg]
            ctgn = ctg.rstrip(args.linename)
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctgn = ctgn[0:9]
            else:
                ctgn = "$" + ctgn + "q"
            img.add(svgwrite.shapes.Rect((xoff+(sc/100),yoff+ypos-3), ((ec-sc)/100,6), stroke='black', stroke_width=1, fill = 'white'))
            if above:
                yt = yoff+ypos-4
                col = "blue" if col == "black" else "black"
            else:
                yt = yoff+ypos+7
            above = not above
            img.add(dwg.text(ctgn, insert=(xoff+(sc/100),yt),fill=col, style="font-size:4"))
            if self.orientation[ctg] == 0:
                direction = ">"
            else:
                direction = "<"
            img.add(dwg.text(direction, insert=(xoff+sc/100,yoff+ypos+2),style="font-size:6"))
        for ctg in sorted(self.contigset_sr, key= lambda x: self.left_coords[x]):
            #print(read)
            sc = self.left_coords[ctg]
            ec = self.right_coords[ctg]
            ctgn = ctg.rstrip(args.linename)
            ctgn = "$" + ctgn + "q"
            img.add(svgwrite.shapes.Rect((xoff+(sc/100),yoff+ypos-3), ((ec-sc)/100,6), stroke='grey', stroke_width=1, fill = 'white'))
            col = "gray"
            if above:
                yt = yoff+ypos-4
            else:
                yt = yoff+ypos+7
            above = not above
            img.add(dwg.text(ctgn, insert=(xoff+(sc/100),yt),fill=col, style="font-size:4"))
            if self.orientation[ctg] == 0:
                direction = ">"
            else:
                direction = "<"
            img.add(dwg.text(direction, insert=(xoff+sc/100,yoff+ypos+2),fill = col, style="font-size:6"))
        return ypos+5

    def set_new_length(self):
        for ctg, coord in self.right_coords.items():
            if self.length < coord:
                self.length = coord

    def move_all_right(self, offset):
        for lr in self.longread_coords:
            self.longread_coords[lr][0] += offset
            self.longread_coords[lr][1] += offset
        for lc in self.left_coords:
            self.left_coords[lc] += offset
        for rc in self.right_coords:
            self.right_coords[rc] += offset

    def add_short_read_contig_left(self, anchor, newctg, distance, orientation):
        try:
            assert(newctg not in self.contigset)
        except AssertionError:
            print("New contig " + str(newctg) + " already exists in " + str(id(self)))
        if distance < 0:
            self.right_coords[newctg] = self.left_coords[anchor]
            self.left_coords[newctg] = self.right_coords[newctg] - contigs[newctg] + distance
        else:
            self.right_coords[newctg] = self.left_coords[anchor] - distance
            self.left_coords[newctg] = self.right_coords[newctg] - contigs[newctg]
        contig2scaffold[newctg] = [id(self)]
        self.orientation[newctg] = orientation
        self.contigset_sr.add(newctg)
        if self.left_coords[newctg] < 0:
            offset = -self.left_coords[newctg]
            self.move_all_right(offset)
        self.set_new_length()

    def add_short_read_contig_right(self, anchor, newctg, distance, orientation):
        try:
            assert(newctg not in self.contigset)
        except AssertionError:
            print("New contig " + str(newctg) + " already exists in " + str(id(self)))
        if distance < 0:
            self.left_coords[newctg] = self.right_coords[anchor]
            self.right_coords[newctg] = self.left_coords[newctg] + allcontigs[newctg] + distance
        else:
            self.left_coords[newctg] = self.right_coords[anchor] + distance
            self.right_coords[newctg] = self.left_coords[newctg] + contigs[newctg]
        contig2scaffold[newctg] = [id(self)]
        self.orientation[newctg] = orientation
        self.contigset_sr.add(newctg)
        self.set_new_length()

    def get_ctg_len(self, ctg):
        if ctg in self.left_coords and ctg in self.right_coords:
            return self.right_coords[ctg] - self.left_coords[ctg]
        else:
            return -1

    def is_on_left_edge(self, ctg):
        if self.left_coords[ctg] < 30:
            return True
        else:
            return False

    def is_on_right_edge(self, ctg):
        if self.length - self.right_coords[ctg] < 30:
            return True
        else:
            return False
    
    def merge_sr(self, ctg1, ctg2, distance):
        try:
            assert(ctg1 in self.contigset or ctg1 in self.contigset_sr)
        except AssertionError:
            print( ctg1 + " is not in " + id(self))
        scaf1 = scaffolds[contig2scaffold[ctg1][0]]
        scaf2 = scaffolds[contig2scaffold[ctg2][0]]

        for rid, read in scaf2.lr_info.items():
            self.lr_info[rid] = read
        for rid, read in scaf2.sr_info.items():
            self.sr_info[rid] = read

        offset = self.right_coords[ctg1] + distance
            
        for ctg,coord in scaf2.left_coords.items():
            self.left_coords[ctg] = coord + offset
        for ctg,coord in scaf2.right_coords.items():
            self.right_coords[ctg] = coord + offset
        for ctg in scaf2.contigset:
            contig2scaffold[ctg] = [id(self)]
            self.contigset.add(ctg)
        for ctg in scaf2.contigset_sr:
            contig2scaffold[ctg] = [id(self)]
            self.contigset_sr.add(ctg)
        for lr,coords in scaf2.longread_coords.items():
            self.longread_coords[lr] = [coords[0] + offset, coords[1] + offset]
        for ctg,ori in scaf2.orientation.items():
            self.orientation[ctg] = ori
        Scaffold.nr_of_scaffolds -= 1
        self.set_new_length()
        del(scaffolds[id(scaf2)])



    def merge(self,scaf2):
        for rid, read in scaf2.lr_info.items():
            self.lr_info[rid] = read
        for rid, read in scaf2.sr_info.items():
            self.sr_info[rid] = read
        same_ctgs = self.contigset.intersection(scaf2.contigset)
        same_orientation= 0
        different_orientation = 0
        for ctg in same_ctgs:
            if self.orientation[ctg] != scaf2.orientation[ctg]:
                different_orientation += 1
            else:
                same_orientation += 1
        if different_orientation > 0 and same_orientation > 0 :
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are oriented differently in the two scaffolds.")
            #print(self.lr_info)
            #print(scaf2.lr_info)
            for ctg in same_ctgs:
                print(ctg + ": " + str(self.orientation[ctg]) + "  " + str(scaf2.orientation[ctg]))
            return    
        if different_orientation > same_orientation : 
            pass
            
        sorted_same_ctgs1 = sorted(same_ctgs, key = lambda item: self.left_coords[item])
        sorted_same_ctgs2 = sorted(same_ctgs, key = lambda item: scaf2.left_coords[item])
        if "-".join(sorted_same_ctgs1) != "-".join(sorted_same_ctgs2):
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are not ordered the same way in the two scaffolds.")
            for ctg in scaf2.contigset:
                try:
                    contig2scaffold[ctg].remove(id(scaf2))
                except ValueError:
                    print("id not found: " + str(id(scaf2)))
                    print("for contig: " + str(ctg))
                    print(contig2scaffold[ctg])
            del(scaffolds[id(scaf2)])
            print(sorted_same_ctgs1)
            print(sorted_same_ctgs2)
            #self.print_contig_sequence()
            #scaf2.print_contig_sequence()
            return
        # This whole thing is not overly complicated but certainly tedious
        # First the scaffold that's more to the left is worked on, this is easier as coordinates don't change much in this scaffold

        lctg = sorted_same_ctgs1[0]
        # smaller means there is less DNA to the left so the scaffold is more right than the other, 
        # length of the contig is important when only part of the contig is mapped 
        if self.is_on_left_edge(lctg) or scaf2.is_on_left_edge(lctg):
            if self.left_coords[lctg] + self.get_ctg_len(lctg)  < scaf2.left_coords[lctg] + scaf2.get_ctg_len(lctg): 
                lscaf = scaf2
                rscaf = self
            else:
                lscaf = self
                rscaf = scaf2
        else:
            if self.left_coords[lctg] < scaf2.left_coords[lctg]:
                lscaf = scaf2
                rscaf = self
            else:
                lscaf = self
                rscaf = scaf2
        nlongread_coords = {}
        for lrid, lrc in lscaf.longread_coords.items():
            nlongread_coords[lrid] = lrc
        lsorted_ctgs= sorted(lscaf.contigset, key = lambda item: lscaf.left_coords[item])
        rsorted_ctgs= sorted(rscaf.contigset, key = lambda item: rscaf.left_coords[item])
        #print("Left scaffold: " + str(id(lscaf)))
        #print("Right scaffold: " + str(id(rscaf)))


        def has_similar_mapping_length(scaf1, ctg1, scaf2, ctg2):
            tolerance = 0.2
            l1 = scaf1.get_ctg_len(ctg1)
            l2 = scaf2.get_ctg_len(ctg2)
            if  l1 == -1:
                return False
            if l2 == -1:
                return False
            if l1 > l2:
                return False if l1 > l2* (1+tolerance) else True
            else:
                return False if l2 > l1* (1+tolerance) else True
        for ctg in same_ctgs:
            lonedge = lsorted_ctgs[0] == ctg or lsorted_ctgs[-1] == ctg
            ronedge = rsorted_ctgs[0] == ctg or rsorted_ctgs[-1] == ctg
            if not (lonedge or ronedge or has_similar_mapping_length(lscaf, ctg, rscaf, ctg)):
                print("WARNING: merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contig " + ctg + " has very different mappging lengths in the two scaffolds.")
        
        last_right_coord = 0
        nleft_coords = {}
        nright_coords = {}
        norientation = {}
        ncontigset = self.contigset.copy()
        last_common_ctg = "nope"

        for ctg in lsorted_ctgs:
            norientation[ctg] = lscaf.orientation[ctg]
            lonedge = lsorted_ctgs[0] == ctg or lsorted_ctgs[-1] == ctg
            #ronedge = rsorted_ctgs[0] == ctg or rsorted_ctgs[-1] == ctg
            #print("lsc: " + lsorted_ctgs[0])
            if lsorted_ctgs[0] == ctg:
                nleft_coords[ctg] = lscaf.left_coords[ctg]
                nright_coords[ctg] = lscaf.right_coords[ctg]
                last_right_coord = lscaf.right_coords[ctg]
                #print(ctg)
                #print(nleft_coords[ctg])
                #print(nright_coords[ctg])
            else: 
                lctgl = lscaf.right_coords[ctg] - lscaf.left_coords[ctg]
                ld = lscaf.left_coords[ctg] - last_right_coord
                #print("lctgl: " + str(lctgl))
                #print("ld: " + str(ld))
                rctgl = lctgl
                rd = ld
                if ctg in same_ctgs:
                    rctgl = rscaf.right_coords[ctg] - rscaf.left_coords[ctg]
                    #rd = rscaf.left_coords[ctg] - last_right_coord
                nleft_coords[ctg] = last_right_coord + ld #(ld if ld < rd else rd)
                nright_coords[ctg] = nleft_coords[ctg] + (lctgl if lctgl > rctgl else rctgl)
                last_right_coords = nright_coords[ctg]
                #print(ctg)
                #print(nleft_coords[ctg])
                #print(nright_coords[ctg])
            #if ctg in same_ctgs:
            #    last_common_ctg = ctg
            #else:
            #    last_common_ctg = "nope"
        first_anchor = "nope"
        for ctg in rsorted_ctgs:
            if ctg in same_ctgs:
                if first_anchor == "nope":
                    first_anchor = ctg
                    if rscaf.is_on_left_edge(ctg):
                        offset = lscaf.right_coords[ctg] - rscaf.right_coords[ctg]
                    else:
                        offset = lscaf.left_coords[ctg] - rscaf.left_coords[ctg]
        for lrid, lrc in rscaf.longread_coords.items():
            nlongread_coords[lrid] = [lrc[0] + offset , lrc[1] + offset]
        last_common_ctg = "nope"
        for ctg in rsorted_ctgs:
            norientation[ctg] = rscaf.orientation[ctg]
            if ctg in same_ctgs:
                last_common_ctg = ctg
                continue
            if last_common_ctg == "nope":
                nright_coords[ctg] = nleft_coords[first_anchor] - (rscaf.left_coords[first_anchor] - rscaf.right_coords[ctg])
                nleft_coords[ctg] = nright_coords[ctg] - rscaf.get_ctg_len( ctg)
            else:
                nleft_coords[ctg] = nright_coords[last_common_ctg] + (rscaf.left_coords[ctg] - rscaf.right_coords[last_common_ctg])
                nright_coords[ctg] = nleft_coords[ctg] + (rscaf.right_coords[ctg] - rscaf.left_coords[ctg])
                
        # collect everything in this scaffold and get rid of scaf2
        llastbit = lscaf.length - lscaf.right_coords[lsorted_ctgs[-1]]
        rlastbit = rscaf.length - rscaf.right_coords[rsorted_ctgs[-1]]
        llength = nright_coords[lsorted_ctgs[-1]] + llastbit
        rlength = nright_coords[rsorted_ctgs[-1]] + rlastbit
        self.longread_coords = nlongread_coords
        self.length = max(llength, rlength)
        self.contigset = self.contigset.union(scaf2.contigset)
        self.left_coords = nleft_coords
        self.right_coords = nright_coords
        self.orientation = norientation
        Scaffold.nr_of_scaffolds -= 1
        for ctg in scaf2.contigset:
            try:
                contig2scaffold[ctg].remove(id(scaf2))
            except ValueError:
                print("id not found: " + str(id(scaf2)))
                print("for contig: " + str(ctg))
                print(contig2scaffold[ctg])
            if id(self) not in contig2scaffold[ctg]:
                contig2scaffold[ctg].append(id(self))
        del(scaffolds[id(scaf2)])
                



    # After merging two scaffolds the coordinate systems have nothing to do with reality anymore
    # The sequences are not taken into account for the merging procedure (the script is called 'less_naive' not 'intricate')
    def merge_contigcoords(self, scaf2):
        for rid, read in scaf2.lr_info.items():
            self.lr_info[rid] = read
        for rid, read in scaf2.sr_info.items():
            self.sr_info[rid] = read
        same_ctgs = self.contigset.intersection(scaf2.contigset)
        #print(same_ctgs)
        same_orientation= 0
        different_orientation = 0
        for ctg in same_ctgs:
            if self.orientation[ctg] != scaf2.orientation[ctg]:
                different_orientation += 1
            else:
                same_orientation += 1
        if different_orientation > 0 and same_orientation > 0 :
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are oriented differentely in the two scaffolds.")
            for ctg in same_ctgs:
                print(ctg + ": " + str(self.orientation[ctg]) + "  " + str(scaf2.orientation[ctg]))
            return    
        if different_orientation > same_orientation : 
            scaf2.turn_around()
            
        offsets = []
        for ctg in same_ctgs:
            offsets.append((scaf2.left_coords[ctg]-self.left_coords[ctg]) - (scaf2.left_coords_contig[ctg] - self.left_coords_contig[ctg]))
            offsets.append((scaf2.right_coords[ctg]-self.right_coords[ctg]) - (scaf2.right_coords_contig[ctg] - self.right_coords_contig[ctg]))
        #print(offsets)
        sorted_same_ctgs1 = sorted(same_ctgs, key = lambda item: self.left_coords[item])
        sorted_same_ctgs2 = sorted(same_ctgs, key = lambda item: scaf2.left_coords[item])
        if "-".join(sorted_same_ctgs1) != "-".join(sorted_same_ctgs2):
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are not ordered the same way in the two scaffolds.")
            self.print_contig_sequence()
            scaf2.print_contig_sequence()
            return
        # This whole thing is not overly complicated but certainly tedious
        # First the scaffold that's more to the left is worked on, this is easier as coordinates don't change much in this scaffold
        lctg = sorted_same_ctgs1[0]
        if self.left_coords[lctg] < scaf2.left_coords[lctg]: # smaller means there is less DNA to the left so the scaffold is more right than the other
            lscaf = scaf2
            rscaf = self
        else:
            lscaf = self
            rscaf = scaf2
        lsorted_ctgs= sorted(lscaf.contigset, key = lambda item: lscaf.left_coords[item])
        rsorted_ctgs= sorted(rscaf.contigset, key = lambda item: rscaf.left_coords[item])

        # The premise for this extension is, the more of a contig can be mappped, the better. 
        # The mapped portion of the contigs is extended, if more of the contig is mapped on the other scaffold
        offset = 0
        nleft_coords = {}
        nright_coords = {}
        nleft_coords_contig = {}
        nright_coords_contig = {}
        norientation = {}
        ncontigset = self.contigset.copy()
        # First the left scaffold is added
        for ctg in lsorted_ctgs:
            norientation[ctg] = lscaf.orientation[ctg]
            if not ctg in same_ctgs:
                nleft_coords[ctg] = lscaf.left_coords[ctg] + offset
                nright_coords[ctg] = lscaf.right_coords[ctg] + offset
                nleft_coords_contig[ctg] = lscaf.left_coords_contig[ctg]
                nright_coords_contig[ctg] = lscaf.right_coords_contig[ctg]
                # potentially new contigs are added to the self-contigset, as self is the scaffold that will ultimately be changed
                ncontigset.add(ctg)
            else: 
                delta_left_coord_contig = lscaf.left_coords_contig[ctg] - rscaf.left_coords_contig[ctg] if lscaf.left_coords_contig[ctg] > rscaf.left_coords_contig[ctg] else 0
                delta_right_coord_contig = rscaf.right_coords_contig[ctg] - lscaf.right_coords_contig[ctg] if lscaf.right_coords_contig[ctg] < rscaf.right_coords_contig[ctg]  else 0
                nleft_coords[ctg] = lscaf.left_coords[ctg]-delta_left_coord_contig + offset
                nright_coords[ctg] = lscaf.right_coords[ctg]-delta_right_coord_contig + offset
                nleft_coords_contig[ctg] = lscaf.left_coords_contig[ctg] - delta_left_coord_contig
                nright_coords_contig[ctg] = lscaf.right_coords_contig[ctg] + delta_right_coord_contig
                offset += delta_left_coord_contig + delta_right_coord_contig

        # For the right scaffold these helper functions are needed
        def get_anchors(scaf1, current_ctg):
            lanchor = None
            ranchor = None
            sorted_ctgs = sorted(scaf1.contigset, key = lambda item: scaf1.left_coords[item])
            ctg_index = sorted_ctgs.index(current_ctg)
            idx = ctg_index - 1
            while idx >= 0: 
                nctg = sorted_ctgs[idx]
                if nctg in same_ctgs:
                    lanchor = nctg
                    break
                idx -= 1
            idx = ctg_index + 1
            while idx < len(sorted_ctgs):
                nctg = sorted_ctgs[idx]
                if nctg in same_ctgs:
                    ranchor = nctg
                    break
                idx += 1
            return (lanchor, ranchor)
        
        def add_with_left_anchor(anchor, ctg, this_scaf):
            delta = this_scaf.left_coords[ctg] - this_scaf.right_coords[anchor] 
            delta += nright_coords_contig[anchor] - this_scaf.right_coords_contig[anchor]
            nleft_coords[ctg] = nright_coords[anchor] + delta
            nright_coords[ctg] = nleft_coords[ctg] + (this_scaf.right_coords_contig[ctg] - this_scaf.left_coords_contig[ctg])
            nleft_coords_contig[ctg] = this_scaf.left_coords_contig[ctg]
            nright_coords_contig[ctg] = this_scaf.right_coords_contig[ctg]

        def add_with_right_anchor(anchor, ctg, this_scaf):
            delta = this_scaf.left_coords[anchor] - this_scaf.right_coords[ctg] 
            delta += nleft_coords_contig[anchor] - this_scaf.left_coords_contig[anchor]
            nright_coords[ctg] = nleft_coords[anchor] - delta
            nleft_coords[ctg] = nright_coords[ctg] - (this_scaf.right_coords_contig[ctg] - this_scaf.left_coords_contig[ctg])
            nleft_coords_contig[ctg] = this_scaf.left_coords_contig[ctg]
            nright_coords_contig[ctg] = this_scaf.right_coords_contig[ctg]

        # Now the right scaffold is added    
        for ctg in rsorted_ctgs:
            if ctg in same_ctgs:
                pass # this has been taken care of in the left scaffold
            else: # contig exclusive to the right scaffold
                ncontigset.add(ctg)
                norientation[ctg] = rscaf.orientation[ctg]
                lanchor, ranchor = get_anchors(rscaf,ctg) 
                if lanchor and ranchor:
                    if rscaf.left_coords[ctg] - rscaf.right_coords[lanchor] < rscaf.left_coords[ranchor] - rscaf.right_coords[ctg]:
                        add_with_left_anchor(lanchor, ctg, rscaf)
                    else:
                        add_with_right_anchor(ranchor, ctg, rscaf)
                elif lanchor:
                    add_with_left_anchor(lanchor, ctg, rscaf)
                else:
                    add_with_right_anchor(ranchor, ctg, rscaf)

        # Calculate new length
        distr = rscaf.length - rscaf.right_coords[rsorted_ctgs[-1]]
        distl = lscaf.length - lscaf.right_coords[lsorted_ctgs[-1]]
        nsorted_ctgs = sorted(ncontigset, key = lambda item: nright_coords[item])
        last_ctg = nsorted_ctgs[-1]
        if last_ctg in lscaf.contigset:
            if last_ctg in rscaf.contigset:
                self.length = nright_coords[last_ctg] + distr if distr > distl else distl
            else:
                self.length = nright_coords[last_ctg] + distl
        else:
            self.length = nright_coords[last_ctg] + distr

        # collect everything in this scaffold and get rid of scaf2
        self.contigset = ncontigset
        self.left_coords = nleft_coords
        self.right_coords = nright_coords
        self.left_coords_contig = nleft_coords_contig
        self.right_coords_contig = nright_coords_contig
        self.orientation = norientation
        Scaffold.nr_of_scaffolds -= 1
        #self.print_contig_sequence()
        #scaf2.print_contig_sequence()
        for ctg in scaf2.contigset:
            contig2scaffold[ctg].remove(id(scaf2))
        del(scaffolds[id(scaf2)])

    # Do not use this!
    # Extracting the sequence is difficult because the sequence data may not be available to you locally.
    # You should write your own method for this. This is my own special case that's why this method has _th as suffix
    # All these paths are hardcoded because the raw data probably won't move around much. 
    def get_sequence_th(self):
        actions = []
        status = "in_lr"
        last_change_pos = 1
        ctgnr = 0
        for pos in range(0,self.length):
            rpos = pos + 1
            if self.coords[pos] != ("","",""):
                (kind, ctg, ctgcoords) = self.coords[pos]
                if kind == "start" and status == "in_lr":
                    status = "in_ctg"
                    last_change_pos = rpos
                    ctgstart = ctgcoords
                elif kind == "end" and status == "in_ctg":
                    status = "in_lr"
                    last_change_pos = rpos
                    action = [ctg,ctgstart,ctgcoords]
                    actions.append(action)
                elif kind == "end-start" and status == "in_ctg":
                    ctg1, ctg2 = ctg.split(",")
                    ctgcoord1, ctgcoord2 = ctgcoords.split(",")
                    last_change_pos = rpos
                    action = [ctg1,ctgstart,ctgcoord1]
                    actions.append(action)
                    ctgstart = ctgcoord2
                elif kind == "start" and status == "in_ctg":
                    status = "overlap"
                    overlap_ctg = ctg
                    overlap_coord = ctgcoords
                    overlap_pos = rpos
                elif kind == "end" and status == "overlap":
                    action = [ctg,ctgstart,ctgcoords]
                    actions.append(action)
                    ctgstart = overlap_coord + rpos - overlap_pos
                    status = "in_ctg"
                elif kind == "end" and status == "in_lr":
                    ctgstart = ctgcoords
                    status = "in_ctg_rev"
                elif kind == "start" and status == "in_ctg_rev":
                    action = [ctg,ctgcoords,ctgstart,"reverse-complementary"]
                    actions.append(action)
                    status = "in_lr"
                else:
                    print("Status \"" + status + "\" and kind \"" + kind + "\" not implemented yet.")
                    print(self.lr_info)
                    self.print_coords()
        print(actions)
        
                

    
    @classmethod
    def init_from_LR(cls,lr):
        newinst = cls()
        newinst.turned_around = False
        newinst.get_sequence_fragments = []
        newinst.scaf_info = []
        newinst.lr_info = dict()
        newinst.lr_info[lr[0]] = lr[1]
        orientation0 = 0
        orientation1 = 0
        newinst.orientation= dict()
        newinst.left_coords = dict()
        newinst.left_coords_contig = dict()
        newinst.right_coords = dict()
        newinst.right_coords_contig = dict()
        newinst.coords = [("","","")]*(lr[1]["length"]+1)
        newinst.length = lr[1]["length"]
        newinst.longread_coords = dict()
        newinst.longread_coords[lr[0]] = [1, lr[1]["length"]]
        for part in lr[1]["maps"]:
            ctg = part["contig"]
            if ctg.endswith(args.linename):
                newinst.contigset.add(ctg)
                if ctg in contigs:
                    del(contigs[ctg])
            # Put information about the mapped contigs in the corresponding data structures
            try: 
                assert(not ctg in newinst.left_coords)
            except AssertionError:
                if ctg.endswith(args.linename):
                    print("Contig " + ctg + " already exists in left_corrds. Probably the contig is in read " + str(lr[0]) + " more than once.")
                    continue
            newinst.left_coords[ctg] = part["scr"]
            newinst.right_coords[ctg] = part["ecr"]
            newinst.orientation[ctg] = part["strand"]

            if ctg.startswith("chr"):
                continue

            if part["strand"] == 1: 
                # workaround for misleading coordinates in erates file
                newinst.left_coords_contig[ctg] = part["lenc"] - part["ecc"]
                newinst.right_coords_contig[ctg] = part["lenc"] - part["scc"]
                startcr = part["ecr"]
                endcr = part["scr"]
            else:
                newinst.left_coords_contig[ctg] = part["scc"]
                newinst.right_coords_contig[ctg] = part["ecc"]
                startcr = part["scr"]
                endcr = part["ecr"]

            if newinst.coords[startcr] != ("","",""):
                if newinst.coords[startcr][0] == "start":
                    print("Problem with coordinates, start")
                elif newinst.coords[startcr][0] == "end":
                    ctgstring = newinst.coords[startcr][1] + "," + ctg
                    ctgcoords = str(newinst.coords[startcr][2]) + "," + str(newinst.right_coords_contig[ctg])
                    newinst.coords[startcr] = ("end-start", ctgstring,ctgcoords)
                else:
                    print("Problem with coordinates, end-start")
            else:
                newinst.coords[startcr] = ("start",ctg,newinst.left_coords_contig[ctg])
            
            if newinst.coords[endcr] != ("","",""):
                if newinst.coords[endcr][0] == "start":
                    ctgstring = ctg + "," + newinst.coords[endcr][1]
                    ctgcoords = str(newinst.right_coords_contig[ctg]) + "," + str(newinst.coords[endcr][2])
                    newinst.coords[endcr] = ("end-start", ctgstring, ctgcoords)
                elif newinst.coords[endcr][0] == "end":
                    print("Problem with coordinates, end")
                else:
                    print("Problem with coordinates, end-start")
            else:
                newinst.coords[endcr] = ("end",ctg, newinst.right_coords_contig[ctg])
            # checking whether the contig is already part of another scaffold happens elsewhere
            if ctg.endswith(args.linename):
                #try:
                #    assert(not ctg in contig2scaffold)
                #except AssertionError:
                #    print("Contig " + ctg + " already exists in contig2scaffold. The reference " + str(contig2scaffold[ctg]) + " will be overwritten.")
                if ctg in contig2scaffold:
                    contig2scaffold[ctg].append(id(newinst))
                else:
                    contig2scaffold[ctg] = [id(newinst)]
            # The orientation of the read is needed.
            # Contigs are somewhat well defined with respect to
            # their orientation, so a majority vote seems appropriate.
            if part["strand"] == 0: 
                orientation0 +=1
                newinst.orientation[ctg] = 0
            else:
                orientation1 +=1
                newinst.orientation[ctg] = 1
        # turn scaffold around
        if orientation0 < orientation1:
            newinst.turn_around()    
        
        
        newinst.idx = id(newinst)
        return newinst
    
    #add_seq_info(self):

# nanopore read info
with open(args.efile) as f:
    for line in f:
        sline = line.split()
        rid = sline[0]
        ctg = sline[1]
        strand = int(sline[8])
        scr = int(sline[5])
        ecr = int(sline[6])
        lenr = int(sline[7])
        scc = int(sline[9])
        ecc = int(sline[10])
        lenc = int(sline[11])
        payload = {"contig":ctg,"strand":strand,"scr":scr,"ecr":ecr,"scc":scc,"ecc":ecc,"lenc":lenc}
        if rid in blacklist:
            if blacklist[rid] == ctg:
                continue
        if rid in reads:
            reads[rid]["maps"].append(payload)
        else:
            reads[rid] = {}
            reads[rid]["length"] = int(sline[7])
            reads[rid]["maps"] = [payload]

scaffolds = {}
# get interesting reads
# and sort contigs by left coordinate
greadst = {}
for rid in reads:
    counter = 0
    for item in reads[rid]["maps"]:
        if item["contig"].endswith(args.linename) and item["contig"] != "1036QBL":
            counter +=1
            if counter >= args.mincontigs:
                greadst[rid] = reads[rid]
                break

for rid in greadst:
    nscaff = Scaffold.init_from_LR((rid,reads[rid]))
    #nscaff.get_sequence_th()
    #nscaff.add_seq_info()
    scaffolds[id(nscaff)] = nscaff


toRemove = set()
for idx,scaf in scaffolds.items():
    #find and get rid of conflicting scaffolds
    if scaf.find_conflicts(10):
        toRemove.add(idx)
    pass
    
for idx in toRemove:
    scaffolds[idx].delete()
    



dwg = svgwrite.Drawing(args.SVG,size=(u'1700', u'2200'), profile='full')
yp = 10
xtext = 10
xpad = 20
dwg.add(dwg.text("10000 bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 10
dwg.add(dwg.line((xpad, yp), ( xpad + 10000/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
for i in range(1,10):
    dwg.add(dwg.line( (xpad + (10000/100)/10 * i, yp-1), (xpad + (10000/100)/10 * i, yp+1), stroke=svgwrite.rgb(0,0,0,'%')))
dwg.add(dwg.line((xpad + 10000/100, yp-2), ( xpad +10000/100, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
yp += 10
dwg.add(dwg.text("100k bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 10
dwg.add(dwg.line((xpad, yp), ( xpad + 100000/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
for i in range(1,10):
    dwg.add(dwg.line( (xpad + (100000/100)/10 * i, yp-1), (xpad + (100000/100)/10 * i, yp+1), stroke=svgwrite.rgb(0,0,0,'%')))
dwg.add(dwg.line((xpad + 100000/100, yp-2), ( xpad +100000/100, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
yp += 10
g1M = dwg.defs.add(dwg.g(id='g003'))
g1M.add(dwg.text("1M bases", insert=( xpad, yp), fill='black', style="font-size:7"))
yp += 10
g1M.add(dwg.line((xpad, yp), ( xpad + 1000000/100, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
g1M.add(dwg.line((xpad, yp-7), ( xpad , yp+7), stroke=svgwrite.rgb(0, 0, 0, '%')))
for i in range(1,10):
    g1M.add(dwg.line( (xpad + (1000000/100)/10 * i, yp-5), (xpad + (1000000/100)/10 * i, yp+5000), stroke=svgwrite.rgb(0,0,0,'%')))
g1M.add(dwg.line((xpad + 1000000/100, yp-7), ( xpad +1000000/100, yp+7), stroke=svgwrite.rgb(0, 0, 0, '%')))
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
        if len(contig2scaffold[contig]) > 1 and len(contig2scaffold[contig]) < 100: # 1036QBL is a problem (139 reads)
            scaf1 = scaffolds[contig2scaffold[contig][0]]
            scaf2 = scaffolds[contig2scaffold[contig][1]]
            scaf1.merge(scaf2)
            break

print("Nr. of scaffolds: " + str(len(scaffolds) + len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")

def is_rightmost(ctg):
    try:
        assert(len(contig2scaffold[ctg]) == 1)
    except AssertionError:
        print("More than one scaffold for contig " + ctg1)
    scaf = scaffolds[contig2scaffold[ctg][0]]
    if ctg == scaf.get_rightmost_contig():
        return True
    return False

def is_leftmost(ctg):
    scaf = scaffolds[contig2scaffold[ctg][0]]
    if ctg == scaf.get_leftmost_contig():
        return True
    return False

if args.summaryfile:
    print("adding short reads ....")
    somethingHappened = True
    while somethingHappened:
        somethingHappened = False
        for ctg1 in srneighs:
            if ctg1 in contig2scaffold and is_rightmost(ctg1):
                scaf1 = scaffolds[contig2scaffold[ctg1][0]]
                for ctg2 in srneighs[ctg1]["right"]:
                    if ctg2[0] in contig2scaffold and is_leftmost(ctg2[0]):
                        scaf2 = scaffolds[contig2scaffold[ctg2[0]][0]]
                        scaf1.merge_sr(ctg1,ctg2[0],ctg2[1])
                        somethingHappened = True
                    elif ctg2[0] in contig2scaffold:
                        pass
                        #print("Problem merging")
                    else:
                        scaf1.add_short_read_contig_right(ctg1,ctg2[0],ctg2[1], 0)
                        del(contigs[ctg2[0]])
                        somethingHappened = True
        for ctg1 in srneighs:
            if ctg1 in contig2scaffold and is_leftmost(ctg1):
                scaf1 = scaffolds[contig2scaffold[ctg1][0]]
                for ctg2 in srneighs[ctg1]["left"]:
                    if ctg2[0] in contig2scaffold: 
                        pass
                    else:
                        scaf1.add_short_read_contig_left(ctg1,ctg2[0],ctg2[1], 0)
                        del(contigs[ctg2[0]])
                        somethingHappened = True

        

# Draw all scaffolds
for scaf in scaffolds.values():
    yp += scaf.to_SVG(dwg, xpad, yp) + 10

dwg.save()
print("Nr. of scaffolds: " + str(len(scaffolds) + len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")

sys.exit(0)
