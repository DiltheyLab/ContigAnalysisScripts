from itertools import combinations, cycle
from collections import defaultdict
import svgwrite
import sys

def shortname(ctgname):
    if "_" in ctgname:
        return ctgname.split("_")[1]
    else:
        return ctgname

def sniff_format(fileh):
    line = fileh.readline()
    rid, nameorlen = line.split()[0:2]
    try:
        int(nameorlen)
        return "paf"
    except:
        return "erate"
        

class Scaffolds:
    lreads = {}
    line = ""
        
    # init just reads in the input file and stores the longreads
    # to get scaffold objects call construct_scaffolds
    def __init__(self, inputfilename, blacklist, line, whitelist = None):
        self.line = line
        with open(inputfilename) as f:
            pafformat = True if sniff_format(f) == "paf" else False
            for line in f:
                if pafformat:
                    [rid, lenr, scr, ecr, strandstring, ctg, lenc, scc, ecc, nr_matches, block_len, quality] = line.split()[0:12]
                    strand = 0 if strandstring == "+" else 1
                else:
                    [rid, ctg, t2, t3, t4, scr, ecr, lenr, strand, scc, ecc, lenc, t12, t13, t14, t15, t16] = line.split()
                    if "strand" == "1": # workaround for misleading coordinates in erates file
                        tmp = int(lenc) - int(ecc)
                        ecc = int(lenc) - int(scc)
                        scc = tmp
                data = {"name":ctg,"strand":int(strand),"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc)}
                if whitelist:
                    if rid not in whitelist:
                        continue
                if blacklist:
                    if rid in blacklist:
                        if "all" in blacklist[rid]:
                            continue
                        elif ctg in blacklist[rid]:
                            continue
                    if shortname(ctg) in blacklist:
                        continue
                if rid in self.lreads:
                    self.lreads[rid]["maps"].append(data)
                    if int(ecr) > self.lreads[rid]["rm_ecr"]:
                        self.lreads[rid]["rm_ecr"] = int(ecr)
                    if int(scr) < self.lreads[rid]["lm_scr"]:
                        self.lreads[rid]["lm_scr"] = int(scr)
                else:
                    self.lreads[rid] = {}
                    self.lreads[rid]["length"] = int(lenr)
                    self.lreads[rid]["maps"] = [data]
                    self.lreads[rid]["rm_ecr"] = int(ecr)
                    self.lreads[rid]["lm_scr"] = int(scr)

    
    def filter_contigcounts(self, nr):
        toremove = set()
        for rid,read in self.lreads.items():
            counter = 0
            for item in read["maps"]:
                if self.line in item["name"]:
                    counter +=1
            
            if counter < nr:
                toremove.add(rid)
        for rid in toremove:
            del(self.lreads[rid])

    def filter_small_contigs(self, size):
        for rid,read in self.lreads.items():
            toremove = []
            for item in read["maps"]:
                if item["ecc"] - item["scc"] < size:
                    toremove.append(item)
            for item in toremove:
                self.lreads[rid]["maps"].remove(item)
    
    def construct_scaffolds(self, contigs):
        scaffolds = {}
        for rid in self.lreads:
            nscaff = Scaffold.init_from_LR((rid,self.lreads[rid]),self.line, contigs )
            scaffolds[id(nscaff)] = nscaff
        return scaffolds

    def sort_contigs_in_reads(self):
        # sort contigs by left coordinate
        for rid in self.lreads:
            self.lreads[rid]["maps"] = sorted(self.lreads[rid]["maps"], key = lambda x: x["scr"])

    def turn_longreads_around(self):
        for rid in self.lreads:
            count0 = 0
            count1 = 0
            for contig in self.lreads[rid]["maps"]:
                if contig["strand"] == 1:
                    count1 += 1
                else:
                    count0 += 1
            if count1 > count0:
                length = self.lreads[rid]["length"]
                for contig in self.lreads[rid]["maps"]:
                    tmp = contig["scr"]
                    contig["scr"] = length - contig["ecr"] 
                    contig["ecr"] = length - tmp
                    contig["strand"] = 0 if contig["strand"] == 1 else 1
                
        

class Scaffold:
    cluster_counter = 0
    linename = None
    in_mergefile = False
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
    name = ""
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
        self.turned_around = False
        #lr_info[lr[0]] = lr[1]
        self.in_mergefile = False


    
    def print_contigset(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: int(item.rstrip(self.linename)))
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
    def turn_around(self, contigs):
        self.turned_around = True
        new_left_coords = dict()
        new_right_coords = dict()
        new_left_coords_contig = dict()
        new_right_coords_contig = dict()
        new_orientation = dict()
        
        for contig in self.contigset:
            new_right_coords[contig] = self.length - self.left_coords[contig] 
            new_left_coords[contig] = self.length - self.right_coords[contig] 
            new_orientation[contig] = 0 if self.orientation[contig] == 1 else 1
            new_left_coords_contig[contig] = contigs[shortname(contig)] - self.right_coords_contig[contig]
            new_right_coords_contig[contig] = contigs[shortname(contig)] - self.left_coords_contig[contig]
        self.left_coords = new_left_coords
        self.right_coords = new_right_coords
        self.left_coords_contig = new_left_coords_contig
        self.right_coords_contig = new_right_coords_contig
        self.orientation = new_orientation
        #self.coords.reverse()
    def remove_revcomp_contigs(self):
        for ctg in self.contigset.copy():
            if self.orientation[ctg] == 1:
                del(self.right_coords[ctg])
                del(self.left_coords[ctg])
                del(self.right_coords_contig[ctg])
                del(self.left_coords_contig[ctg])
                del(self.orientation[ctg])
                self.contigset.remove(ctg)

    def get_contained_contigs(self):
        output = set()
        for ctg1,ctg2 in combinations(self.contigset,2):
            if self.left_coords[ctg2] < self.left_coords[ctg1] and self.right_coords[ctg2] > self.right_coords[ctg1]:
                output.add((ctg1,ctg2))
            elif self.left_coords[ctg1] < self.left_coords[ctg2] and self.right_coords[ctg1] > self.right_coords[ctg2]:
                output.add((ctg2,ctg1))
        return output


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

    # length that takes into account that the end of the scaffold could be extended by a mapping contig
    def get_virtual_length(self, contigs):
        rctg = self.get_rightmost_contig()
        rest =  contigs[rctg] - self.right_coords_contig[rctg]
        return max(self.length ,self.right_coords[rctg] + rest)

    def find_conflicts(self):
        sortedcontigs = sorted(self.contigset, key = lambda item: self.left_coords[item])
        sortedcontigs2 = sorted(self.contigset, key = lambda item: self.right_coords[item])
        scstring1 = "-".join(sortedcontigs)
        scstring2 = "-".join(sortedcontigs2)
        if scstring1 != scstring2:
            print("Conflict! \n" + scstring1 + "\n" + scstring2 + " . Removing scaffold.")
            return True
        return False

    def remove_overlapped_contigs(self, contigs):
        # first remove the overlapped contigs
        leftc = sorted(self.contigset, key = lambda item: self.left_coords[item])
        rightc = sorted(self.contigset, key = lambda item: self.right_coords[item])
        toremove = set()
        for ctg1 in self.contigset:
            for ctg2 in self.contigset:
                if self.left_coords[ctg2] < self.left_coords[ctg1] and self.right_coords[ctg2] > self.right_coords[ctg1]:
                    toremove.add(ctg1)
        for ctg in toremove:
            self.contigset.remove(ctg)
            del(self.left_coords[ctg])
            del(self.right_coords[ctg])
            del(self.left_coords_contig[ctg])
            del(self.right_coords_contig[ctg])
            #contig2scaffold[ctg].remove(id(self))
            print("Removed contig " + str(ctg) + " from " + str(self.name))
        return toremove

    def remove_short_contigs(self, contigs):
        # Now lets get rid of contigs with no space and less than 0.5 mapped length
        toremove = set()
        sctgs = sorted(self.contigset, key= lambda x: self.left_coords[x])
        for ctgl, ctg, ctgr in zip(sctgs[:-2],sctgs[1:-1],sctgs[2:]):
            #print("-".join([ctgl,ctg,ctgr]))
            space_left = self.left_coords[ctg] - self.right_coords[ctgl]
            space_right = self.left_coords[ctgr] - self.right_coords[ctg]
            space = space_left + space_right
            #print(self.right_coords_contig[ctg])
            if self.right_coords_contig[ctg] - self.left_coords_contig[ctg] + space < 0.5 * contigs[shortname(ctg)]:
                toremove.add(ctg)
        for ctg, ctgr in zip(sctgs[:-1],sctgs[1:]):
            space_right = self.left_coords[ctgr] - self.right_coords[ctg]
            rest_right = contigs[shortname(ctg)] - self.right_coords_contig[ctg]
            if rest_right - space_right > 1000:
                toremove.add(ctg)
        for ctgl, ctg in zip(sctgs[:-1],sctgs[1:]):
            space_left = self.left_coords[ctg] - self.right_coords[ctgl]
            rest_left = self.left_coords_contig[ctg]
            if rest_left - space_left > 1000:
                toremove.add(ctg)

        for ctg in toremove:
            self.contigset.remove(ctg)
            del(self.left_coords[ctg])
            del(self.right_coords[ctg])
            del(self.left_coords_contig[ctg])
            del(self.right_coords_contig[ctg])
            #contig2scaffold[ctg].remove(id(self))
            #print("Removed contig " + str(ctg) + " from " + str(self.name) + " due to lack of space.")
        return toremove
        

    # Returns the space in y that it needs (depends on the number of longreads that were merged into this scaffold)
    def to_SVG(self, img, contigs, xoff, yoff,show_lr_ids):
        ypad = 7
        col = "black"
        nr_longreads = len(self.longread_coords)
        ypos = 0
        if show_lr_ids:
            y_space_per_longread = 4
            for lrid, lrc in self.longread_coords.items():
                ylen = nr_longreads * y_space_per_longread - ypos
                rect = img.add(svgwrite.shapes.Rect((xoff+(lrc[0]/100),yoff+ypos), ((lrc[1]-lrc[0])/100,ylen+ypad), stroke='green', stroke_width=1 ))
                rect.fill(color="none").dasharray([2, 2])
                img.add(img.text(lrid, insert=(xoff+(lrc[0]/100),yoff+ypos-1),fill="green", style="font-size:2"))
                ypos += y_space_per_longread
        ypos += ypad
            
        img.add(svgwrite.shapes.Line((xoff, yoff+ypos), ( xoff + self.length/100, yoff+ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))

        ctg_y_drawsize = 8
        ctg_y_halfdrawsize = ctg_y_drawsize/2
        ctg_relative_positions = cycle([-ctg_y_halfdrawsize-1, ctg_y_halfdrawsize+3, -ctg_y_halfdrawsize-4, ctg_y_halfdrawsize+6])

        gradient_idc = 0
        for ctg in sorted(self.contigset, key= lambda x: self.left_coords[x]):
            #print(read)
            sc = self.left_coords[ctg]
            ec = self.right_coords[ctg]
            scc = self.left_coords_contig[ctg]
            ecc = self.right_coords_contig[ctg]
            ctgn = ctg.rstrip(self.linename)
            ctgn2 = ctg.rstrip("rc").rstrip(self.linename)
            #ctg = read[0]
            if ctg.startswith("chr"):
                ctgn = ctgn[0:9]
            else:
                ctgn = "$" + ctgn + "q"

            gradient_idc += 1
            gradient_id = self.name + "_" + str(gradient_idc)
            lineargrad = img.defs.add(svgwrite.gradients.LinearGradient(id=gradient_id , x1=-scc/(ecc-scc), x2=1+(contigs[shortname(ctg)]-ecc)/(ecc-scc), y1=0, y2=0))
            if ctgn2.endswith("0") or ctgn2.endswith("1"):
                col1 = "#FF0000"
                col2 = "#0000FF"
            elif ctgn2.endswith("2") or ctgn2.endswith("3"):
                col1 = "#FFE119"
                col2 = "#000000"
            elif ctgn2.endswith("4") or ctgn2.endswith("5"):
                col1 = "#911eb4"
                col2 = "#a9a9a9"
            elif ctgn2.endswith("6") or ctgn2.endswith("7"):
                col1 = "#000075"
                col2 = "#f58231"
            elif ctgn2.endswith("8") or ctgn2.endswith("9"):
                col1 = "#e6beff"
                col2 = "#808000"
                
            lineargrad.add_stop_color("0%",col1)
            lineargrad.add_stop_color("50%","#FFFFFF")
            lineargrad.add_stop_color("100%",col2)

            #leftclip = scc/100
            #rightclip = (contigs[shortname(ctg)]-ecc)/100
            img.add(svgwrite.shapes.Rect((xoff+(sc/100),yoff+ypos-ctg_y_halfdrawsize), ((ec-sc)/100,ctg_y_drawsize), stroke='black', stroke_width=1, fill = 'url(#'+gradient_id + ')'))
            #g.add(svgwrite.shapes.Rect((xpad+((xoffset+sc)/100),ypad+ypos-6), ((ec-sc)/100,12), stroke='black', stroke_width=1, fill='url(#'+str(gradient_idc)+')'))
            yt = yoff + ypos + next(ctg_relative_positions)
            img.add(img.text(ctgn, insert=(xoff+(sc/100),yt),fill=col, style="font-size:3"))
            if self.orientation[ctg] == 0:
                direction = ">"
            else:
                direction = "<"
            img.add(img.text(direction, insert=(xoff+sc/100,yoff+ypos+2),style="font-size:6"))
        for ctg in sorted(self.contigset_sr, key= lambda x: self.left_coords[x]):
            #print(read)
            sc = self.left_coords[ctg]
            ec = self.right_coords[ctg]
            ctgn = ctg.rstrip(self.linename)
            ctgn = "$" + ctgn + "q"
            img.add(svgwrite.shapes.Rect((xoff+(sc/100),yoff+ypos-ctg_y_halfdrawsize), ((ec-sc)/100,ctg_y_drawsize), stroke='grey', stroke_width=1, fill = 'white'))
            col = "gray"
            yt = yoff + ypos + next(ctg_relative_positions)
            img.add(img.text(ctgn, insert=(xoff+(sc/100),yt),fill=col, style="font-size:3"))
            if self.orientation[ctg] == 0:
                direction = ">"
            else:
                direction = "<"
            img.add(img.text(direction, insert=(xoff+sc/100,yoff+ypos+2),fill = col, style="font-size:6"))
        return ypos+7

    # use only when sr_contigs are added
    # long reads can be longer even after contigs are added
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

    def add_short_read_contig_right(self, anchor, newctg, distance, orientation, contigs, mergefile=None):
        try:
            assert(newctg not in self.contigset)
        except AssertionError:
            print("New contig " + str(newctg) + " already exists in " + str(id(self)))
        self.left_coords[newctg] = self.right_coords[anchor] + (contigs[anchor] - self.right_coords_contig[anchor]) + distance 
        self.right_coords[newctg] = self.left_coords[newctg] + contigs[newctg]
        self.left_coords_contig[newctg] = 1
        self.right_coords_contig[newctg] = contigs[newctg]
        #contig2scaffold[newctg] = [id(self)]
        self.orientation[newctg] = orientation
        self.contigset_sr.add(newctg)
        oldlength = self.length
        self.set_new_length()
        newname = "cluster_" + str(Scaffold.cluster_counter) 
        Scaffold.cluster_counter += 1
        #self.left_coords[newctg] = self.right_coords[anchor] + contigs[anchor] - self.right_coords_contig[anchor] + distance
        full_distance = self.right_coords[anchor] + contigs[anchor] - self.right_coords_contig[anchor] + distance 
        if anchor == "2406APD":
            print(self.left_coords[anchor])
            print(self.right_coords[anchor])
            print(self.left_coords[newctg])
            print(self.right_coords[newctg])

        # add merge info to mergefile 
        if mergefile:
            mode = "extension_ctg" 
            with open(mergefile, "a+") as mergef:
                mergef.write("\t".join([mode ,self.name, str(oldlength), str(self.turned_around), newctg, str(contigs[newctg]), str(False), str(full_distance), newname]))
                mergef.write("\n")
        self.name = newname

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
    
    def merge_sr(self, ctg1, ctg2, scaf1, scaf2, distance, contigs, mergefile = None):
        try:
            assert(ctg1 in self.contigset or ctg1 in self.contigset_sr)
        except AssertionError:
            print( ctg1 + " is not in " + id(self))
        #scaf1 = scaffolds[contig2scaffold[ctg1][0]]
        #scaf2 = scaffolds[contig2scaffold[ctg2][0]]
        if id(scaf1) != id(self):
            print("Problem with ids. self :" + str(id(self)) + " scaf1: " + str(id(scaf1)))
        vlength = self.get_virtual_length(contigs)
        e = contigs[ctg1] - scaf1.right_coords_contig[ctg1] + distance + scaf2.left_coords_contig[ctg2]
        offset = scaf1.right_coords[ctg1] + e
        lr_dist = offset - scaf2.left_coords[ctg2]
        nlength = lr_dist + scaf2.length

        for rid, read in scaf2.lr_info.items():
            self.lr_info[rid] = read
        for rid, read in scaf2.sr_info.items():
            self.sr_info[rid] = read

        for ctg,coord in scaf2.left_coords.items():
            self.left_coords[ctg] = coord + offset
        for ctg,coord in scaf2.right_coords.items():
            self.right_coords[ctg] = coord + offset
        for ctg,coord in scaf2.left_coords_contig.items():
            self.left_coords_contig[ctg] = coord 
        for ctg,coord in scaf2.right_coords_contig.items():
            self.right_coords_contig[ctg] = coord 
        #for ctg in scaf2.contigset:
        #    contig2scaffold[ctg] = [id(self)]
        #    self.contigset.add(ctg)
        #for ctg in scaf2.contigset_sr:
        #    contig2scaffold[ctg] = [id(self)]
        #    self.contigset_sr.add(ctg)
        for lr,coords in scaf2.longread_coords.items():
            self.longread_coords[lr] = [coords[0] + offset, coords[1] + offset]
        for ctg,ori in scaf2.orientation.items():
            self.orientation[ctg] = ori
        Scaffold.nr_of_scaffolds -= 1
        #self.set_new_length()
        newname = "cluster_" + str(Scaffold.cluster_counter)
        Scaffold.cluster_counter += 1
        if mergefile:
            mode = "merging" 
            with open(mergefile, "a+") as mergef:
                mergef.write("\t".join([mode ,self.name, str(self.length), str(self.turned_around), scaf2.name, str(scaf2.length), str(scaf2.turned_around), str(lr_dist), newname]))
                mergef.write("\n")
        self.length = nlength
        self.name = newname
        #del(contig2scaffold[ctg2][0])
        #del(scaffolds[id(scaf2)])
        return True


    def get_distance_pairs(self, ctgs):
        dist_matrix = {}
        for ctg1, ctg2  in combinations(ctgs, 2):
            if self.left_coords[ctg1] > self.left_coords[ctg2]:
                lctg, rctg = [ctg1, ctg2]
            else:
                rctg, lctg = [ctg1, ctg2]
            dist = self.left_coords[rctg] - self.left_coords_contig[rctg] - self.left_coords[lctg] + self.left_coords_contig[lctg]
            dist_matrix[(lctg, rctg)] = dist
            dist_matrix[(rctg, lctg)] = -dist
        return dist_matrix

    def merge(self,scaf2, mergefile = None):
        for rid, read in scaf2.lr_info.items():
            self.lr_info[rid] = read
        for rid, read in scaf2.sr_info.items():
            self.sr_info[rid] = read
        same_ctgs = self.contigset & scaf2.contigset
        same_orientation= 0
        different_orientation = 0
        for ctg in same_ctgs:
            if self.orientation[ctg] != scaf2.orientation[ctg]:
                different_orientation += 1
            else:
                same_orientation += 1
        if different_orientation > 0 and same_orientation > 0 :
            print("Problem merging " + str(self.name) + " and " + str(scaf2.name) + ". Contigs are oriented differently in the two scaffolds.")
            for ctg in same_ctgs:
                print(ctg + ": " + str(self.orientation[ctg]) + "  " + str(scaf2.orientation[ctg]))
            return False
        if different_orientation > same_orientation : 
            pass

        # Sanity Check (hard): contig ordering -> contig distances
        def get_pairs_with_proper_distances(scaf1, scaf2, ctgs, tolerance):
            pairs1 = scaf1.get_distance_pairs(ctgs)
            #print(pairs1)
            pairs2 = scaf2.get_distance_pairs(ctgs)
            #print(pairs2)
            pairs = set()
            for pair in pairs1.keys():
                ldist = abs(pairs1[pair]) if abs(pairs1[pair]) > abs(pairs1[pair]) else abs(pairs2[pair])
                maxdist = max(int(ldist*tolerance), 4000)
                if pairs1[pair] > pairs2[pair]:
                    if pairs1[pair] < pairs2[pair] + maxdist:
                        pairs.add(pair)
                else:
                    if pairs2[pair] < pairs1[pair] + maxdist:
                        pairs.add(pair)
            return pairs
            
        #if not same_ctgs:
        #    print(self.name + " " + scaf2.name + " have no common contigs")
        #print(same_ctgs)
        goodpairs = get_pairs_with_proper_distances(self, scaf2, same_ctgs, 0.6)
        problempairs = set(combinations(same_ctgs,2)) - goodpairs
        if problempairs:
        #    pass
            print("problempairs: " + str(problempairs))
            #print("goodpairs: " + str(goodpairs))
            #print(self.get_distance_pairs(same_ctgs))
            #print(scaf2.get_distance_pairs(same_ctgs))
            print(self.contigset & scaf2.contigset)
            print(self.left_coords)
            print(self.left_coords_contig)
            print(scaf2.left_coords)
            print(scaf2.left_coords_contig)
            for ctg1, ctg2 in problempairs:
                print(self.get_distance_pairs(set([ctg1,ctg2])))
                print(scaf2.get_distance_pairs(set([ctg1,ctg2])))
                #print("\t".join([ctg1, str(self.left_coords[ctg1]),str(self.left_coords_contig[ctg1]), ctg2, str(self.left_coords[ctg2]), str(self.left_coords_contig[ctg2])]))
                #print("\t".join([ctg1, str(scaf2.left_coords[ctg1]),str(scaf2.left_coords_contig[ctg1]), ctg2, str(scaf2.left_coords[ctg2]), str(scaf2.left_coords_contig[ctg2])]))
            print(self.name + "  " + scaf2.name)
            #print(str(self.lr_info.keys()) + "  " + str(scaf2.lr_info.keys()))
            sys.exit()

        sorted_same_ctgs1 = sorted(same_ctgs, key = lambda item: self.left_coords[item])
        sorted_same_ctgs2 = sorted(same_ctgs, key = lambda item: scaf2.left_coords[item])
        if "-".join(sorted_same_ctgs1) != "-".join(sorted_same_ctgs2):
            print("Problem merging " + str(self.name) + " and " + str(scaf2.name) + ". Contigs are not ordered the same way in the two scaffolds.")
            print(sorted_same_ctgs1)
            print(sorted_same_ctgs2)
            self.print_contig_sequence()
            scaf2.print_contig_sequence()
            return False

        # This whole thing is not overly complicated but certainly tedious
        # First the scaffold that's more to the left is worked on, this is easier as coordinates don't change much in this scaffold
        # So the left scaffold is needed which will be determined by the leftmost anchor
        lctg = sorted_same_ctgs1[0]
        if self.right_coords[lctg] + self.right_coords_contig[lctg] < scaf2.right_coords[lctg] + scaf2.right_coords_contig[lctg]:
            rscaf = self
            lscaf = scaf2
        else:
            rscaf = scaf2
            lscaf = self
        
        # smaller means there is less DNA to the left so the scaffold is more right than the other, 
        # length of the contig is important when only part of the contig is mapped 
        # This may be redundant and may have to get cleaned up later on.
        # The distance of the two longreads is needed.
        # As best guess the distance is taken that is given by the common contig
        # that is rightmost on both scaffolds
        rctg = sorted_same_ctgs1[-1]
        distance = (lscaf.left_coords[rctg] - lscaf.left_coords_contig[rctg]) - (rscaf.left_coords[rctg] - rscaf.left_coords_contig[rctg])

        # add merge info to mergefile 
        mode = "incorporation" if distance + rscaf.length < lscaf.length else "extension"
        if mergefile:
            self.in_mergefile = True 
            with open(mergefile, "a+") as mergef:
                mergef.write("\t".join([mode ,lscaf.name, str(lscaf.length), str(lscaf.turned_around), rscaf.name, str(rscaf.length), str(rscaf.turned_around), str(distance), "cluster_" + str(Scaffold.cluster_counter)]))
                mergef.write("\n")

        nlongread_coords = {}
        for lrid, lrc in lscaf.longread_coords.items():
            nlongread_coords[lrid] = lrc
        lsorted_ctgs= sorted(lscaf.contigset, key = lambda item: lscaf.left_coords[item])
        rsorted_ctgs= sorted(rscaf.contigset, key = lambda item: rscaf.left_coords[item])
        

        def has_similar_mapping_length(scaf1, ctg1, scaf2, ctg2):
            tolerance = 0.5
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

        # Sanity Check (soft): mapping lengths
        for ctg in same_ctgs:
            lonedge = lsorted_ctgs[0] == ctg or lsorted_ctgs[-1] == ctg
            ronedge = rsorted_ctgs[0] == ctg or rsorted_ctgs[-1] == ctg
            if not (lonedge or ronedge or has_similar_mapping_length(lscaf, ctg, rscaf, ctg)):
                pass
                #print("WARNING: merging " + str(self.idx) + " and " + str(scaf2.idx) + ". Contig " + ctg + " has very different mapping lengths (" + str(lscaf.get_ctg_len(ctg)) + "/" + str(rscaf.get_ctg_len(ctg)) + ") in the two scaffolds.")
        nleft_coords = {}
        nleft_coords_contig = {}
        nright_coords = {}
        nright_coords_contig = {}
        norientation = {}
        ncontigset = self.contigset.copy()
        first_anchor = None

        for ctg in rscaf.contigset - lscaf.contigset:
            nright_coords_contig[ctg] = rscaf.right_coords_contig[ctg]
            nleft_coords_contig[ctg] = rscaf.left_coords_contig[ctg]
        for ctg in lscaf.contigset:
            nright_coords_contig[ctg] = lscaf.right_coords_contig[ctg]
            nleft_coords_contig[ctg] = lscaf.left_coords_contig[ctg]

        # take same coordinates for contigs on left scaffold
        for ctg in lsorted_ctgs:
            norientation[ctg] = lscaf.orientation[ctg]
            nleft_coords[ctg] = lscaf.left_coords[ctg]
            nright_coords[ctg] = lscaf.right_coords[ctg]
        # take offset for first common contig
        for ctg in rsorted_ctgs:
            if ctg in same_ctgs:
                first_anchor = ctg
                offset = (lscaf.left_coords[ctg] - lscaf.left_coords_contig[ctg]) -  (rscaf.left_coords[ctg] - rscaf.left_coords_contig[ctg])
                break
        # utilize offset for long read info
        for lrid, lrc in rscaf.longread_coords.items():
            nlongread_coords[lrid] = [lrc[0] + offset , lrc[1] + offset]
        last_common_ctg = None
        for ctg in rsorted_ctgs:
            norientation[ctg] = rscaf.orientation[ctg]
            if ctg in same_ctgs:
                last_common_ctg = ctg
                offset = (lscaf.left_coords[ctg] - lscaf.left_coords_contig[ctg]) -  (rscaf.left_coords[ctg] - rscaf.left_coords_contig[ctg])
            else:
                nright_coords[ctg] = rscaf.right_coords[ctg] + offset
                nleft_coords[ctg] = rscaf.left_coords[ctg] + offset

        # Sanity Check (soft): contig spaces
        new_overlapping = 0
        example_problem = None
        con1 = self.get_contained_contigs()
        con2 = scaf2.get_contained_contigs()
                
        # collect everything in this scaffold and get rid of scaf2
        llastbit = lscaf.length - lscaf.right_coords[lsorted_ctgs[-1]]
        rlastbit = rscaf.length - rscaf.right_coords[rsorted_ctgs[-1]]
        llength = nright_coords[lsorted_ctgs[-1]] + llastbit
        rlength = nright_coords[rsorted_ctgs[-1]] + rlastbit
        self.longread_coords = nlongread_coords
        self.length = max(llength, rlength)
        self.contigset = self.contigset.union(scaf2.contigset)
        self.left_coords = nleft_coords
        self.left_coords_contig = nleft_coords_contig
        self.right_coords_contig = nright_coords_contig
        self.right_coords = nright_coords
        self.orientation = norientation
        Scaffold.nr_of_scaffolds -= 1

        connew = self.get_contained_contigs()
        if len(connew) > len(con1 | con2):
            print("Problem merging " + str(self.name) + " and " + str(scaf2.name) + ". Merge has introduced " + str(len(connew)) + " new overlaps.")
            print("Merge by " + str(same_ctgs))
            print(connew - (con1 | con2))

        for ctg in scaf2.contigset:
            if ctg in same_ctgs:
                #self.lef_coords_contig[ctg] = 
                #self.left_coords_contig[ctg] = lscaf.left_coords_contig[ctg] if lscaf.left_coords_contig[ctg] < rscaf.left_coords_contig[ctg] else rscaf.left_coords_contig[ctg]
                #self.right_coords_contig[ctg] = lscaf.right_coords_contig[ctg] if lscaf.right_coords_contig[ctg] > rscaf.right_coords_contig[ctg] else rscaf.right_coords_contig[ctg]
                pass # let's keep the info from the self scaffold
            else:
                self.left_coords_contig[ctg] = scaf2.left_coords_contig[ctg]
                self.right_coords_contig[ctg] = scaf2.right_coords_contig[ctg]
            #if id(self) not in contig2scaffold[ctg]:
            #    contig2scaffold[ctg].append(id(self))
            #if id(scaf2) in contig2scaffold[ctg]:
            #    contig2scaffold[ctg].remove(id(scaf2))
        self.name = "cluster_" + str(Scaffold.cluster_counter) 
        Scaffold.cluster_counter += 1
        #del(scaffolds[id(scaf2)])
        return True

                

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

    # Deprecated: the information is now saved during merging
    # Do not use this!
    # Extracting the sequence is difficult because the sequence data may not be available to you locally.
    # You should write your own method for this. This is my own special case that's why this method has _th as suffix
    # All these paths are hardcoded because the raw data probably won't move around much. 
    def get_sequence_th(self):
        actions = []
        status = "in_lr"
        last_change_pos = 1
        ctgs = []
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
    def init_from_LR(cls,lr,linename,contigs):
        newinst = cls()
        newinst.turned_around = False
        newinst.in_mergefile = False
        newinst.linename = linename
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
            ctg = part["name"]
            if ctg.endswith(newinst.linename):
                if ctg in newinst.contigset:
                    lold = abs(newinst.right_coords[ctg] - newinst.left_coords[ctg])
                    lnew = abs(part["ecr"]-part["scr"])
                    if lnew < lold: # ignore this contig-map, because it is smaller than the previous one
                        continue
                    if lnew < 300:
                        continue
                else:
                    newinst.contigset.add(ctg)
            # Put information about the mapped contigs in the corresponding data structures
            newinst.left_coords[ctg] = part["scr"]
            newinst.right_coords[ctg] = part["ecr"]
            newinst.orientation[ctg] = part["strand"]

            if ctg.startswith("chr"):
                continue

            # checking whether the contig is already part of another scaffold happens during merging procedure
            newinst.left_coords_contig[ctg] = part["scc"]
            newinst.right_coords_contig[ctg] = part["ecc"]
            startcr = part["scr"]
            endcr = part["ecr"]

            # The orientation of the read is needed.
            # Contigs are somewhat well defined with respect to
            # their orientation, so a majority vote seems appropriate.
            if part["strand"] == 0: 
                orientation0 +=1
                newinst.orientation[ctg] = 0
            else:
                orientation1 +=1
                newinst.orientation[ctg] = 1
        # turn scaffold around if needed -> now done on the long read level
        #if orientation0 < orientation1:
        #    newinst.turn_around(contigs)    
        # get rid of wrong contigs
        newinst.remove_revcomp_contigs()
        newinst.name = lr[0]
        newinst.idx = id(newinst)
        return newinst
    #add_seq_info(self):
