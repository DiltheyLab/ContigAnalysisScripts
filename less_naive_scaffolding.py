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
    
    def add_lr_info(self, lr):
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
        print("Turned " + str(id(self)) + " around.")
        new_left_coords = {}
        new_right_coords = {}
        new_orientation = {}
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
        
    
    @classmethod
    def init_from_LR(cls,lr):
        newinst = cls()
        newinst.lr_info[lr[0]] = lr[1]
        orientation0 = 0
        orientation1 = 0
        newinst.left_coords = dict()
        newinst.right_coords = dict()
        newinst.coords = [("","")]*lr[1]["length"]
        newinst.length = lr[1]["length"]
        for part in lr[1]["maps"]:
            ctg = part["contig"]
            if ctg.endswith("QBL"):
                newinst.contigset.add(ctg)
            # the orientation of the read is needed.
            # Contigs are somewhat well defined with respect to
            # their orientation, so a majority vote seems appropriate.
            # This information is used later to turn the scaffold around
            if part["strand"] == 0: 
                orientation0 +=1
                newinst.orientation[ctg] = 0
            else:
                orientation1 +=1
                newinst.orientation[ctg] = 1
            # Put information about the mapped contigs in the corresponding data structures
            try: 
                assert(not ctg in newinst.left_coords)
            except AssertionError:
                if ctg.startswith("QBL"):
                    print("Contig " + ctg + " already exists in left_corrds. Probably the contig is in read " + str(lr[0]) + " more than once.")
            newinst.left_coords[ctg] = part["scr"]
            try:
                assert(not ctg in newinst.right_coords)
            except AssertionError:
                if ctg.startswith("QBL"):
                    print("Contig " + ctg + " already exists in right_corrds. Probably the contig is in read " + str(lr[0]) + " more than once.")
            newinst.right_coords[ctg] = part["ecr"]
            newinst.coords[part["scr"]-1] = ("start",ctg)
            newinst.coords[part["ecr"]-1] = ("end",ctg)
            # checking whether the contig is already part of another scaffold happens elsewhere
            if ctg.startswith("QBL"):
                try:
                    assert(not ctg in contig2scaffold)
                except AssertionError:
                    print("Contig " + ctg + " already exists in contig2scaffold. The reference " + str(contig2scaffold[ctg]) + " will be overwritten.")
                contig2scaffold[ctg] = id(newinst)
        # turn scaffold around
        if orientation0 < orientation1 and orientation1 > 1:
            newinst.turn_around()    
           

        return newinst

    @staticmethod
    def return_contig_set():
        return contigset

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
    #print(reads[rid]["overlaps"])
    nscaff = Scaffold.init_from_LR((rid,reads[rid]))
    #nscaff.add_lr_info((rid,greadst[rid]))
    scaffolds[id(nscaff)] = nscaff
    
    
    
    #soverlaps = sorted(reads[rid]["maps"], key = itemgetter("scr"))
    #greads[rid] = greadst[rid]
    #greads[rid]["maps"]=soverlaps

for idx,scaf in scaffolds.items():
    #print(idx)
    #scaf.print_coords()
    pass

sys.exit(0)


# cluster np-reads 
print("scaffolding long reads ....")
creads = {}
clusternr = 0
while len(greads) > 0:
    clusternr += 1
    current_cluster = {}
    current_contigs = set()
    # take a random read and build a cluster from it
    cr = greads.popitem()
    current_cluster[cr[0]] = cr[1]
    olen = 0
    while len(current_cluster) != olen:
        olen = len(current_cluster)
        for contig in cr[1]["overlaps"]:
            if not contig["contig"].startswith("chr"):
                contigs.pop(contig["contig"], None)
                current_contigs.add(contig["contig"])
                contig2cluster[contig["contig"]] = clusternr
        for readid,readval in greads.items():
            contig_found = False
            for contig in readval["overlaps"]:
                if not contig["contig"].startswith("chr"):
                    if contig["contig"] in current_contigs:
                        contig_found = True
                        current_cluster[readid] = readval
                        cr = (readid, greads.pop(readid))
                        break
                   
            if contig_found:
                break
    creads[clusternr] = current_cluster
print("Nr. of scaffolds: " + str(clusternr+len(contigs)) + " (" + str(clusternr) + " cluster + " + str(len(contigs))+ " contigs)")

# simpler data structure to collect contigs into scaffolds
#scaffolds={}
for i, cluster in creads.items():
    current_contigs = set([])
    for readid,read in cluster.items():
        #print(read)
        for contig in read["overlaps"]:
            if not contig["contig"].startswith("chr"):
                current_contigs.add(contig["contig"])
    scaffolds[i] = current_contigs
#print(scaffolds)

def addcontig(ctg, cluster):
    contig2cluster[ctg] = cluster
    scaffolds[cluster].add(ctg)
    contigs.pop(ctg, None)
    
def mergecluster(cluster1, cluster2):
    for contig in scaffolds[cluster2]:
        contig2cluster[contig] = cluster1
        scaffolds[cluster1].add(contig)
    scaffolds.pop(cluster2)

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
