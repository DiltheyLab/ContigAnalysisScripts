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
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")
parser.add_argument("SVG", help="Scaffolds are drawn to this SVG file")
#parser.add_argument("--maxdev", help="Maximal deviation", type=float, default=2.0)
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=10)

args = parser.parse_args()

reads = {}
greads = {}
contigs = {}
contig2scaffold = {}


for read in SeqIO.parse(args.contigfile, "fasta"):
    contigs[read.id] = len(read.seq)

print("Nr. of scaffolds: " + str(len(contigs)))

class Scaffold:
    sr_info = dict()
    lr_info = dict()
    longread_coords = dict()
    coords = []
    left_coords = dict()
    right_coords = dict()
    orientation = dict()
    contigset = set()
    length = 0
    sequence = ""
    nr_of_scaffolds = 0
    turned_around = False
    idx = ""
    #id = 
    # All contigs will have scaffold coordinates.
    # Before scaffolds are merged all scaffold coordinates 
    # are equal to long read coordinates of the orginal long read 
    # that initialized the scaffold
    
    def __init__(self):
        Scaffold.nr_of_scaffolds += 1
        self.contigset = set()
        #self.turned_around = False
        #lr_info[lr[0]] = lr[1]

    def add_sr_info(self, sr):
        pass
    

    def print_contigset(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: int(item.rstrip("QBL")))
        for contig in sortedcontigs:
            print(contig)
        
    def print_coords(self):
        coordstr = ""
        for pos in self.coords:
            if pos == ("",""):
                coordstr += "."
            else:
                coordstr += "|"
        print(coordstr)

    def print_contig_sequence(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: self.left_coords[item])
        print("-".join(sortedcontigs))
        
    def print_id(self):
        print(self.idx)

    def del_lr_info(self, lrid):
        pass
        # can be complicated
        # probably it's easiest to create a new object
        # with all the reads except for the specified lrid

    def del_sr_info(self, srid):
        pass

    # this signifies a hard turnaround
    # i.e. the coordinate system of the scaffold is reversed
    def turn_around(self):
        #self.turned_around = not self.turned_around
        #print("Turned " + str(id(self)) + " around.")
        new_left_coords = dict()
        new_right_coords = dict()
        new_orientation = dict()
        #if "54QBL" in self.contigset:
        #    print(self.left_coords)
        #    print(self.right_coords)
        for contig in self.contigset:
            new_right_coords[contig] = self.length - self.left_coords[contig] 
            new_left_coords[contig] = self.length - self.right_coords[contig] 
            new_orientation[contig] = 0 if self.orientation[contig] == 1 else 1
        self.left_coords = new_left_coords
        self.right_coords = new_right_coords
        self.orientation = new_orientation
        #if "54QBL" in self.contigset:
        #    print(self.left_coords)
        #    print(self.right_coords)

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

    def get_leftmost(self,ctgset):
        sortedcontigs = sorted(ctgset, key = lambda item: self.left_coords[item])
        return(sortedcontigs[0])


    def find_conflicts(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: self.left_coords[item])
        sortedcontigs2 = sorted(self.contigset, key = lambda item: self.right_coords[item])
        scstring1 = "-".join(sortedcontigs)
        scstring2 = "-".join(sortedcontigs2)
        if scstring1 != scstring2:
            print("Conflict!")
            print(scstring1)
            print(scstring2)

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
            ctgn = ctg.rstrip("QBL")
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctgn = ctgn[0:9]
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
        return ypos+5

    def add_short_read_contig(self, anchor, newctg, distance, orientation):
        pass
        
        


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
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are oriented differentely in the two scaffolds.")
            for ctg in same_ctgs:
                print(ctg + ": " + str(self.orientation[ctg]) + "  " + str(scaf2.orientation[ctg]))
            return    
        if different_orientation > same_orientation : 
            scaf2.turn_around()
            
        sorted_same_ctgs1 = sorted(same_ctgs, key = lambda item: self.left_coords[item])
        sorted_same_ctgs2 = sorted(same_ctgs, key = lambda item: scaf2.left_coords[item])
        if "-".join(sorted_same_ctgs1) != "-".join(sorted_same_ctgs2):
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are not ordered the same way in the two scaffolds.")
            self.print_contig_sequence()
            scaf2.print_contig_sequence()
            return
        # This whole thing is not overly complicated but certainly tedious
        # First the scaffold that's more to the left is worked on, this is easier as coordinates don't change much in this scaffold
        def get_ctg_len(scaf, ctg):
            if ctg in scaf.left_coords and ctg in scaf.right_coords:
                return scaf.right_coords[ctg] - scaf.left_coords[ctg]
            else:
                return -1

        def is_on_left_edge(scaf, ctg):
            if scaf.left_coords[ctg] < 30:
                return True
            else:
                return False
        def is_on_right_edge(scaf, ctg):
            if scaf.length - scaf.right_coords[ctg] < 30:
                return True
            else:
                return False

        lctg = sorted_same_ctgs1[0]
        # smaller means there is less DNA to the left so the scaffold is more right than the other, 
        # length of the contig is important when only part of the contig is mapped 
        if is_on_left_edge(self,lctg) or is_on_left_edge(scaf2,lctg):
            if self.left_coords[lctg] + get_ctg_len(self, lctg)  < scaf2.left_coords[lctg] + get_ctg_len(scaf2, lctg): 
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
            l1 = get_ctg_len(scaf1,ctg1)
            l2 = get_ctg_len(scaf2,ctg2)
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
                if is_on_left_edge(rscaf, ctg):
                    offset = lscaf.right_coords[ctg] - rscaf.right_coords[ctg]
                else:
                    offset = lscaf.left_coords[ctg] - rscaf.left_coords[ctg]
        for lrid, lrc in rscaf.longread_coords.items():
            nlongread_coords[lrid] = (lrc[0] + offset , lrc[1] + offset)
        last_common_ctg = "nope"
        for ctg in rsorted_ctgs:
            norientation[ctg] = rscaf.orientation[ctg]
            if ctg in same_ctgs:
                last_common_ctg = ctg
                continue
            if last_common_ctg == "nope":
                nright_coords[ctg] = nleft_coords[first_anchor] - (rscaf.left_coords[first_anchor] - rscaf.right_coords[ctg])
                nleft_coords[ctg] = nright_coords[ctg] - get_ctg_len(rscaf, ctg)
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
        #if "54QBL" in scaf2.contigset:
        #    print("ok?")
        #    print(scaf2.lr_info.keys())
        #    print(scaf2.left_coords)
        #    print(scaf2.right_coords)
        #if "54QBL" in self.contigset:
        #    print("hwhat?")
        #    print(self.lr_info.keys())
        #    print(self.left_coords)
        #    print(self.right_coords)
        #print("merging " + str(id(self)) + " and " + str(id(scaf2)))
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
                #if "54QBL" in self.contigset:
                #    print("nleft_coords:" + str(nleft_coords))
                #    print("offset:" + str(offset))
                #    print("dlcc: " + str(delta_left_coord_contig))
                #    print("drcc: " + str(delta_right_coord_contig))
                #    print("llcc: " + str(lscaf.left_coords_contig[ctg]))
                #    print("lrcc: " + str(lscaf.right_coords_contig[ctg]))
                #    print("rlcc: " + str(rscaf.left_coords_contig[ctg]))
                #    print("rrcc: " + str(rscaf.right_coords_contig[ctg]))

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
        
    
    @classmethod
    def init_from_LR(cls,lr):
        newinst = cls()
        newinst.lr_info = dict()
        newinst.lr_info[lr[0]] = lr[1]
        orientation0 = 0
        orientation1 = 0
        newinst.orientation= dict()
        newinst.left_coords = dict()
        newinst.left_coords_contig = dict()
        newinst.right_coords = dict()
        newinst.right_coords_contig = dict()
        newinst.coords = [("","")]*lr[1]["length"]
        newinst.length = lr[1]["length"]
        newinst.longread_coords = dict()
        newinst.longread_coords[lr[0]] = (1, lr[1]["length"])
        for part in lr[1]["maps"]:
            ctg = part["contig"]
            if ctg.endswith("QBL"):
                newinst.contigset.add(ctg)
                if ctg in contigs:
                    del(contigs[ctg])
            # Put information about the mapped contigs in the corresponding data structures
            try: 
                assert(not ctg in newinst.left_coords)
            except AssertionError:
                if ctg.endswith("QBL"):
                    print("Contig " + ctg + " already exists in left_corrds. Probably the contig is in read " + str(lr[0]) + " more than once.")
                    continue
            newinst.left_coords[ctg] = part["scr"]
            newinst.left_coords_contig[ctg] = part["scc"]
            newinst.right_coords[ctg] = part["ecr"]
            newinst.right_coords_contig[ctg] = part["ecc"]
            newinst.coords[part["scr"]-1] = ("start",ctg)
            newinst.coords[part["ecr"]-1] = ("end",ctg)
            # checking whether the contig is already part of another scaffold happens elsewhere
            if ctg.endswith("QBL"):
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
        if item["contig"].endswith("QBL") and item["contig"] != "1036QBL":
            counter +=1
            if counter == 2:
                greadst[rid] = reads[rid]
                break

for rid in greadst:
    nscaff = Scaffold.init_from_LR((rid,reads[rid]))
    scaffolds[id(nscaff)] = nscaff

for idx,scaf in scaffolds.items():
    scaf.find_conflicts()



dwg = svgwrite.Drawing(args.SVG,size=(u'1700', u'4600'), profile='full')
yp = 10
xtext = 10
xpad = 20
dwg.add(dwg.text("10000 bases", insert=( xpad, yp), fill='black', style="font-size:7"))
dwg.add(dwg.line((xpad, yp+4), ( xpad + 10000/100, yp+4), stroke=svgwrite.rgb(0, 0, 0, '%')))
dwg.add(dwg.line((xpad, yp+2), ( xpad , yp+6), stroke=svgwrite.rgb(0, 0, 0, '%')))
for i in range(1,10):
    dwg.add(dwg.line( (xpad + (10000/100)/10 * i, yp+3), (xpad + (10000/100)/10 * i, yp+5), stroke=svgwrite.rgb(0,0,0,'%')))
dwg.add(dwg.line((xpad + 10000/100, yp+2), ( xpad +10000/100, yp+6), stroke=svgwrite.rgb(0, 0, 0, '%')))
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
            #yp += 20
            #dwg.add(dwg.text(id(scaf1), insert=(xtext, yp+1), fill='black', style="font-size:7"))
            #scaf1.to_SVG(dwg, xpad, yp)
            #yp += 20
            #dwg.add(dwg.text(id(scaf2), insert=(xtext, yp+1), fill='black', style="font-size:7"))
            #scaf2.to_SVG(dwg, xpad, yp)
            #print(contig2scaffold[contig])
            scaf1.merge(scaf2)
            #yp += 20
            #dwg.add(dwg.text(id(scaf1), insert=(xtext, yp+1), fill='black', style="font-size:7"))
            #scaf1.to_SVG(dwg, xpad, yp)
            #yp += 20
            break
print("Nr. of scaffolds: " + str(len(scaffolds) + len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")

print("adding short reads ....")
with open(args.summaryfile) as f:
    for line in f:
        sline = line.split()
        ctg1 = sline[0].split("_")[0].strip("+").strip("-")
        ctg2 = sline[0].split("_")[1].strip("+").strip("-")
        if sline[1] == "NA":
            continue
        if int(sline[2]) < args.mindepth:
            continue
        #moddist = float(sline[1])
        if ctg1 in contig2scaffold:
            if ctg2 in contig2scaffold:
                if contig2scaffold[ctg1] != contig2scaffold[ctg2]:
                    pass # TODO implement short read merging
                else
                    pass
            else:
                scaf = scaffolds[contig2scaffold[ctg1]]
                scaf.add_shortread_contig(ctg1, ctg2, )
            #print("cluster of ctg1 (" + ctg1 + "): " + str(contig2cluster[ctg1]))
        elif ctg2 in contig2cluster:
            addcontig(ctg1, contig2cluster[ctg2])
        else:
            clusternr += 1
            #print("new cluster: " + str(clusternr))
            scaffolds[clusternr] = set([ctg1, ctg2])
            contig2cluster[ctg1] = clusternr
            contig2cluster[ctg2] = clusternr
            contigs.pop(ctg1, None)
            contigs.pop(ctg2, None)

# Draw all scaffolds
for scaf in scaffolds.values():
    yp += scaf.to_SVG(dwg, xpad, yp) + 10

dwg.save()

sys.exit(0)

