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
parser.add_argument("--mindepth", help="Minimal depth", type=int, default=20)

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
        if "54QBL" in self.contigset:
            print(self.left_coords)
            print(self.right_coords)
        for contig in self.contigset:
            new_right_coords[contig] = self.length - self.left_coords[contig] 
            new_left_coords[contig] = self.length - self.right_coords[contig] 
            new_orientation[contig] = 0 if self.orientation[contig] == 1 else 1
        self.left_coords = new_left_coords
        self.right_coords = new_right_coords
        self.orientation = new_orientation
        if "54QBL" in self.contigset:
            print(self.left_coords)
            print(self.right_coords)

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

    def to_SVG(self, img, xoff, yoff):
        if "54QBL" in self.contigset:
            print(self.left_coords)
            print(self.right_coords)
        ypad = 10
        img.add(dwg.line((xoff, yoff+ypad), ( xoff + self.length/100, yoff+ypad), stroke=svgwrite.rgb(0, 0, 0, '%')))
        above = True
        col = "black"
        for ctg in self.contigset:
            #print(read)
            sc = self.left_coords[ctg]
            ec = self.right_coords[ctg]
            ctgn = ctg.rstrip("QBL")
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctgn = ctgn[0:9]
            img.add(svgwrite.shapes.Rect((xoff+(sc/100),yoff+ypad-3), ((ec-sc)/100,6), stroke='black', stroke_width=1, fill = 'white'))
            if above:
                yt = yoff+ypad-4
                col = "blue" if col == "black" else "black"
            else:
                yt = yoff+ypad+7
            above = not above
            img.add(dwg.text(ctgn, insert=(xoff+sc/100,yt),fill=col, style="font-size:4"))
            if self.orientation[ctg] == 0:
                direction = ">"
            else:
                direction = "<"
            img.add(dwg.text(direction, insert=(xoff+sc/100,yoff+ypad+2),style="font-size:6"))

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
        lctg = sorted_same_ctgs1[0]
        if self.left_coords[lctg] < scaf2.left_coords[lctg]: # smaller means there is less DNA to the left so the scaffold is more right than the other
            lscaf = scaf2
            rscaf = self
        else:
            lscaf = self
            rscaf = scaf2
        lsorted_ctgs= sorted(lscaf.contigset, key = lambda item: lscaf.left_coords[item])
        rsorted_ctgs= sorted(rscaf.contigset, key = lambda item: rscaf.left_coords[item])

        def get_ctg_len(self, ctg):
            if ctg in self.left_coords and ctg in self.right_coords:
                return self.right_coords[ctg] - self.left_coords[ctg]
            else:
                return -1

        def has_similar_mapping_length(self, ctg1, scaf2, ctg2):
            tolerance = 0.2
            l1 = self.get_ctg_len(ctg1)
            l2 = scaf2.get_ctg_len(ctg2)
            if  l1 == -1:
                return False
            if l2 == -1:
                return False
            if l1 > l2:
                if l1 > l2* (1+tolerance):
                    return False
                else:
                    return True
            else:
                if l2 > l1* (1+tolerance):
                    return False
                else:
                    return True

            
            
        
        # First the left coordinate in the LR of the right scaffold is needed
        # therefore the left-most contig that has the same length (with tolerance) in the scaffolds is needed
        
        # If the right scaffold does not extend further right than the left scaffold
        # we only check whether the contig coordinates are similar enough,
        # and if any contigs or on one scaffold but not the other
        if rsorted_ctgs[-1] in lsorted_ctgs:



    # After merging two scaffolds the coordinate systems have nothing to do with reality anymore
    # The sequences are not taken into account for the merging procedure (the script is called 'less_naive' not 'intricate')
    def merge_contigcoords(self, scaf2):
        if "54QBL" in scaf2.contigset:
            print("ok?")
            print(scaf2.lr_info.keys())
            print(scaf2.left_coords)
            print(scaf2.right_coords)
        if "54QBL" in self.contigset:
            print("hwhat?")
            print(self.lr_info.keys())
            print(self.left_coords)
            print(self.right_coords)
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
                if "54QBL" in self.contigset:
                    print("nleft_coords:" + str(nleft_coords))
                    print("offset:" + str(offset))
                    print("dlcc: " + str(delta_left_coord_contig))
                    print("drcc: " + str(delta_right_coord_contig))
                    print("llcc: " + str(lscaf.left_coords_contig[ctg]))
                    print("lrcc: " + str(lscaf.right_coords_contig[ctg]))
                    print("rlcc: " + str(rscaf.left_coords_contig[ctg]))
                    print("rrcc: " + str(rscaf.right_coords_contig[ctg]))

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
        newinst.longread_coords[id(newinst)] = (1, lr[1]["length"])
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
            greadst[rid] = reads[rid]
            break

for rid in greadst:
    nscaff = Scaffold.init_from_LR((rid,reads[rid]))
    scaffolds[id(nscaff)] = nscaff

for idx,scaf in scaffolds.items():
    scaf.find_conflicts()




# cluster np-reads 
print("scaffolding long reads ....")
creads = {}
clusternr = 0
olen_scaf = len(scaffolds)+1
while len(scaffolds)-olen_scaf != 0:
    olen_scaf = len(scaffolds)
    for contig in contig2scaffold:
        if len(contig2scaffold[contig]) > 1 and len(contig2scaffold[contig]) < 100: # 1036QBL is a problem (139 reads)
            #print(contig2scaffold[contig])
            scaffolds[contig2scaffold[contig][0]].merge(scaffolds[contig2scaffold[contig][1]])
            print("Nr. of scaffolds: " + str(len(scaffolds) + len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")
            break


# Draw SVG
ypos = 0
xtext = 10
xpad = 200
ypad = 10
dwg = svgwrite.Drawing(args.SVG,size=(u'1700', u'4600'), profile='full')
csize = 0
toggle = False
for scaf in scaffolds.values():
    #print(scaf.orientation)
    dwg.add(dwg.text(id(scaf), insert=(xtext, ypad+ypos+1), fill='black', style="font-size:7"))
    ypos += 20
    scaf.to_SVG(dwg, xpad, ypos)

dwg.save()

sys.exit(0)

