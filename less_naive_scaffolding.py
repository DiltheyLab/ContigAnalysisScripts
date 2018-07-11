from argparse import ArgumentParser
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from itertools import combinations


parser = ArgumentParser()
parser.add_argument("efile", help="Error rate file")
parser.add_argument("summaryfile", help="Contig Distance Summary file")
parser.add_argument("contigfile", help="Contig File")
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
    coords = []
    left_coords = dict()
    right_coords = dict()
    orientation = dict()
    contigset = set()
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
        for contig in self.contigset:
            new_right_coords[contig] = self.length - self.left_coords[contig] 
            new_left_coords[contig] = self.length - self.right_coords[contig] 
            new_orientation[contig] = 0 if self.orientation[contig] == 1 else 1
        self.left_coords = new_left_coords
        self.right_coords = new_right_coords
        self.orientation = new_orientation

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

    # After merging two scaffolds the coordinate systems have nothing to do with reality anymore
    # The sequences are not taken into account for the merging procedure (the script is called 'less_naive' not 'intricate')
    def merge(self, scaf2):
        for rid, read in scaf2.lr_info.items():
            self.lr_info[rid] = read
        for rid, read in scaf2.sr_info.items():
            self.sr_info[rid] = read
        same_orientation = 0
        same_ctgs = self.contigset.intersection(scaf2.contigset)
        #print(same_ctgs)
        same_orientation=0
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
        sortedctgs1 = sorted(same_ctgs, key = lambda item: self.left_coords[item])
        sortedctgs2 = sorted(same_ctgs, key = lambda item: scaf2.left_coords[item])
        if "-".join(sortedctgs1) != "-".join(sortedctgs2):
            print("Problem merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contigs are not ordered the same way in the two scaffolds.")
            return
        
        for nr,ctg in sortedctgs1.items():
        # the easier case 
            if (scaf2.left_coords[ctg]-self.left_coords[ctg]) - (scaf2.left_coords_contig[ctg] - self.left_coords_contig[ctg]) > 0: 
                if scaf2.left_coords_contig[ctg] > 
        
        #print(self.get_leftmost(same_ctgs))
        #print(scaf2.get_leftmost(same_ctgs))
        
        
        
    
    @classmethod
    def init_from_LR(cls,lr):
        newinst = cls()
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
        for part in lr[1]["maps"]:
            ctg = part["contig"]
            if ctg.endswith("QBL"):
                newinst.contigset.add(ctg)
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
        if orientation0 < orientation1 and orientation1 > 1:
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
        if item["contig"].endswith("QBL"):
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
print("Nr. of scaffolds: " + str(clusternr+len(contigs)) + " (" + str(clusternr) + " cluster + " + str(len(contigs))+ " contigs)")

sys.exit(0)



# very lenient clustering of short reads
print("scaffolding short reads ....")
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
        if ctg1 in contig2cluster:
            if ctg2 in contig2cluster:
                if contig2cluster[ctg1] != contig2cluster[ctg2]:
                    #print("merging clusters " + str(contig2cluster[ctg1]) + " and " + str(contig2cluster[ctg2]))
                    mergecluster(contig2cluster[ctg1], contig2cluster[ctg2])
            else:
                addcontig(ctg2, contig2cluster[ctg1])
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
            
            
print("Nr. of scaffolds: " + str(len(scaffolds)+len(contigs)) + " (" + str(len(scaffolds)) + " cluster + " + str(len(contigs))+ " contigs)")
#print(contigs)
for scaf in scaffolds:
    print(scaffolds[scaf])
    sys.exit(0)
def get_rightmost_contig(scaf_id):
    print(creads[scaf_id])
    if len(creads[scaf_id]) > 1:
        print("not implemented yet")
    else:
        rid, read = creads[scaf_id].items()
        for ov in read["overlaps"]:
            print(ov)

#get_rightmost_contig(1)
#sys.exit(0)
