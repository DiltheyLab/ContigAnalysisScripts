from itertools import combinations, cycle, product
from collections import defaultdict, Counter, deque, OrderedDict
from statistics import mean
import svgwrite
from svgwrite.container import Group
import sys
import re

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
        
complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
def revcomp(instring):
    outstring = ""
    for char in instring:
        outstring += complement[char]
    return outstring[::-1]

class Longreads(object):
    lreads = {}
    cellline = ""
    ctg2lreads = defaultdict(set)
    contig_lengths = {}

    def delete(self, rid):
        self.lreads.pop(rid)
        for ctg, lrids in self.ctg2lreads.items():
            if rid in lrids:
                lrids.remove(rid)

    def __str__(self):
        out = ""
        for lread in self.lreads.values():
            #out += str(lread)
            out += "Length: " + str(lread["length"]) + "\n"
            for mapper in lread["maps"]:
                out += "\t" + str(mapper) + "\n"
        return out

    def print_ids(self):
        for rid in self.lreads.keys():
            print(rid)
            #print(lread["name"])
        
    # init just reads in the input file and stores the longreads
    # to get scaffold objects call construct_scaffolds
    def __init__(self, inputfiles, blacklist, line, whitelist_lreads = None):
        self.cellline = line
        for inputfilename in inputfiles:
            with open(inputfilename) as f:
                pafformat = True if sniff_format(f) == "paf" else False
            with open(inputfilename) as f:
                for line in f:
                    if pafformat:
                        [rid, lenr, scr, ecr, strandstring, ctg, lenc, scc, ecc, nr_matches, block_len, quality] = line.split()[0:12]
                        strand = 0 if strandstring == "+" else 1
                        eratequality = 100
                    else:
                        [rid, ctg, t2, eratequality, t4, scr, ecr, lenr, strandstring, scc, ecc, lenc, t12, t13, t14, t15, t16] = line.split()
                        strand = 0 if strandstring == "0" else 1
                        if strand == 1: # workaround for misleading coordinates in erates file
                            tmp = int(lenc) - int(ecc)
                            ecc = int(lenc) - int(scc)
                            scc = tmp
                    if self.cellline not in ctg:
                        continue
                    self.contig_lengths[ctg] = int(lenc)
                    data = {"name":ctg,"strand":strand,"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc), "quality": float(eratequality)/100.0}
                    # Filtering 
                    if whitelist_lreads:
                        if rid not in whitelist_lreads:
                            continue
                    if blacklist:
                        if rid in blacklist:
                            if "all" in blacklist[rid]:
                                continue
                            elif ctg in blacklist[rid]:
                                continue
                        if shortname(ctg) in blacklist:
                            if (int(ecc)-int(scc))/int(lenc) < blacklist[ctg]:
                                continue
                    self.ctg2lreads[ctg].add(rid)
                    if rid in self.lreads:
                        self.lreads[rid]["ctgset"].add(ctg)
                        self.lreads[rid]["mapsc"][ctg].append(data)
                        self.lreads[rid]["maps"].append(data)
                        if int(ecr) > self.lreads[rid]["rm_ecr"]:
                            self.lreads[rid]["rm_ecr"] = int(ecr)
                        if int(scr) < self.lreads[rid]["lm_scr"]:
                            self.lreads[rid]["lm_scr"] = int(scr)
                    else:
                        self.lreads[rid] = {}
                        self.lreads[rid]["ctgset"] = set()
                        self.lreads[rid]["ctgset"].add(ctg)
                        self.lreads[rid]["length"] = int(lenr)
                        # maps is just a list. mapsc allows for easy access to all contigs with a certain name
                        self.lreads[rid]["maps"] = [data]
                        self.lreads[rid]["mapsc"] = defaultdict(list)
                        self.lreads[rid]["mapsc"][ctg].append(data)
                        self.lreads[rid]["rm_ecr"] = int(ecr)
                        self.lreads[rid]["lm_scr"] = int(scr)
                        self.lreads[rid]["fname"] = inputfilename.split("/")[-1]

    @classmethod
    def init_from_reverse_paf(cls, inputfilename):
        newinst = cls([],None,"n.a.")
        with open(inputfilename) as f:
            for line in f:
                [ctg, lenc, scc, ecc, strandstring, tmp1, lenr, scr, ecr, nr_matches, block_len, quality] = line.split()[0:12]
                m = re.search('[A-Z]+', ctg)
                rid = m.group(0)
                strand = 0 if strandstring == "+" else 1
                data = {"name":ctg,"strand":strand,"scr":int(scr),"ecr":int(ecr),"scc":int(scc),"ecc":int(ecc),"lenc":int(lenc), "quality": 1.0}
                if rid in newinst.lreads:
                    newinst.lreads[rid]["ctgset"].add(ctg)
                    newinst.lreads[rid]["mapsc"][ctg].append(data)
                    newinst.lreads[rid]["maps"].append(data)
                    if int(ecr) > newinst.lreads[rid]["rm_ecr"]:
                        newinst.lreads[rid]["rm_ecr"] = int(ecr)
                    if int(scr) < newinst.lreads[rid]["lm_scr"]:
                        newinst.lreads[rid]["lm_scr"] = int(scr)
                else:
                    newinst.lreads[rid] = {}
                    newinst.lreads[rid]["ctgset"] = set()
                    newinst.lreads[rid]["ctgset"].add(ctg)
                    newinst.lreads[rid]["length"] = int(lenr)
                    # maps is just a list. mapsc allows for easy access to all contigs with a certain name
                    newinst.lreads[rid]["maps"] = [data]
                    newinst.lreads[rid]["mapsc"] = defaultdict(list)
                    newinst.lreads[rid]["mapsc"][ctg].append(data)
                    newinst.lreads[rid]["rm_ecr"] = int(ecr)
                    newinst.lreads[rid]["lm_scr"] = int(scr)
                    newinst.lreads[rid]["fname"] = inputfilename.split("/")[-1]
        return newinst

    @classmethod
    def init_from_dict(cls, lrdict, line, contig_lengths, lrids = None):
        newinst = cls([],None,line)
        newinst.cellline = line
        newinst.lrids = lrids
        newinst.contig_lengths = contig_lengths
        newinst.lreads = defaultdict(dict)
        for rid, longread in lrdict.items():
            newinst.lreads[rid]["rm_ecr"] = None
            newinst.lreads[rid]["lm_scr"] = None
            newinst.lreads[rid]["ctgset"] = set()
            newinst.lreads[rid]["maps"] = []
            newinst.lreads[rid]["mapsc"] = defaultdict(list)
            for contig in lrdict[rid]["maps"]:
                ctgn = contig["name"]
                newinst.ctg2lreads[ctgn].add(rid)
                newinst.lreads[rid]["ctgset"].add(ctgn)
                newinst.lreads[rid]["mapsc"][ctgn].append(contig)
                newinst.lreads[rid]["maps"].append(contig)
                if (not newinst.lreads[rid]["lm_scr"]) or contig["scr"] < newinst.lreads[rid]["lm_scr"]:
                    newinst.lreads[rid]["lm_scr"] = contig["scr"]
                if (not newinst.lreads[rid]["rm_ecr"]) or contig["ecr"] > newinst.lreads[rid]["rm_ecr"]:
                    newinst.lreads[rid]["rm_ecr"] = contig["ecr"]
            newinst.lreads[rid]["length"] = newinst.lreads[rid]["rm_ecr"]
            
        return newinst

    def filter_whitelist_ctgs(self, whitelist_ctgs):
        for rid in list(self.lreads.keys()):
            if self.lreads[rid]["ctgset"] & whitelist_ctgs:
                continue
            del(self.lreads[rid])

    
    def filter_contigcounts(self, nr):
        toremove = set()
        for rid,read in self.lreads.items():
            counter = 0
            for item in read["maps"]:
                if self.cellline in item["name"]:
                    counter +=1
            if counter < nr:
                toremove.add(rid)
        for rid in toremove:
            for item in self.lreads[rid]["maps"]:
                if rid in self.ctg2lreads[item["name"]]:
                    self.ctg2lreads[item["name"]].remove(rid)
            del(self.lreads[rid])

    def remove_contig_from_read(self, rid, ctg):
        self.lreads[rid]["maps"].remove(ctg)
        ctgn = ctg["name"]
        self.lreads[rid]["mapsc"][ctgn].remove(ctg)
        if not self.lreads[rid]["mapsc"][ctgn]: # if there is no contig with this name left in the read
            self.ctg2lreads[ctgn].remove(rid)
            self.lreads[rid]["ctgset"].remove(ctgn)

    def filter_low_quality_contigs(self, quality):
        for rid, lread in self.lreads.items():
            for ctg in lread["maps"]:
                if ctg["quality"] < quality:
                    self.remove_contig_from_read(rid, ctg)

    def get_distance_matrix(self, matrixfilename, srdata = None):
        with open(matrixfilename, "w+") as out:
            for rid,read in self.lreads.items():
                for idx1, ctg1 in enumerate(read["maps"][:-1]):
                    for ctg2 in read["maps"][idx1+1:]:
                        dist = ctg2["scr"]-ctg2["scc"] - (ctg1["ecr"] + ctg1["lenc"] - ctg1["ecc"])
                        out.write("\t".join([ctg1["name"], ctg2["name"], str(dist)])+ "\n" )
            if srdata:
                for ctgs, dist in srdata.items():
                    out.write(ctgs[0] + "\t" + ctgs[1] + "\t" + str(dist) + "\n")
                

    def filter_contigs_by_coverage(self, coverage, ignore_ends=False, verbose=False):
        for rid, lread in self.lreads.items():
            for ctg in lread["maps"]:
                if (ctg == lread["maps"][0] or ctg == lread["maps"][-1]) and ignore_ends:
                    continue
                if (ctg["ecc"] - ctg["scc"])/ctg["lenc"] < coverage:
                    self.remove_contig_from_read(rid, ctg)
                    if verbose:
                        print("Removed " + str(ctg["name"]) + " from " + str(rid))


    # assume contigs are sorted
    def filter_overlapped_contigs(self, verbose=False):
        for rid,read in self.lreads.items():
            toremove = []
            last_end = 0
            for ctg in sorted(read["maps"], key = lambda x: x["scr"]):
                if ctg["ecr"] < last_end:
                    toremove.append(ctg)
                else:
                    last_end = ctg["ecr"]
            for item in toremove:
                self.remove_contig_from_read(rid, item)
                if verbose:
                    print("Removed " + item["name"] + " from " + rid + ".")

    def filter_small_contigs(self, size):
        for rid,read in self.lreads.items():
            toremove = []
            for item in read["maps"]:
                if item["ecc"] - item["scc"] < size:
                    toremove.append(item)
            for item in toremove:
                self.remove_contig_from_read(rid, item)

    def filter_small_double_contigs(self, ctglengths, fraction, verbose=False):
        contigcounts = Counter()
        self.sort_by_starts()
        for rid,read in self.lreads.items():
            for ctg in read["maps"]:
                contigcounts[ctg["name"]] += 1
        multis = set([x for x in contigcounts.keys() if contigcounts[x] > 1])
        for rid,read in self.lreads.items():
            for ctg in read["maps"][1:-1]: # first and last contig of each read get a pass
                if ctg["name"] in multis:
                    if (ctg["ecc"] - ctg["scc"])/ctglengths[ctg["name"]] < fraction:
                        self.remove_contig_from_read(rid, ctg)
                        contigcounts[ctg["name"]] -= 1
                        if verbose:
                            pass
                            #print("Removed " + ctg["name"] + " from " + str(rid) + "\t" + str(contigcounts[ctg["name"]]) + " remaining.")
        if verbose:
            for ctgn in set([x for x in contigcounts.keys() if contigcounts[x] != 1]):
                print(ctgn + ": " + str(contigcounts[ctgn]))

    def filter_reverse_double_contigs(self, ctglengths, verbose=False):
        contigcounts = Counter()
        for rid,read in self.lreads.items():
            for ctg in read["maps"]:
                contigcounts[ctg["name"]] += 1
        multis = set([x for x in contigcounts.keys() if contigcounts[x] > 1])
        for rid,read in self.lreads.items():
            for ctg in read["maps"]:
                if ctg["name"] in multis:
                    if ctg["strand"] == 1:
                        self.remove_contig_from_read(rid, ctg)
                        contigcounts[ctg["name"]] -= 1
        if verbose:
            #known_ctgs = set([x+"APD" for x in ["401","411","388", "339","58","1497","2328", "2185", "558", "419", "1697", "52", "47", "162", "1004", "370","462","407"]])
            #print(known_ctgs)
            for ctgn in set([x for x in contigcounts.keys() if contigcounts[x] != 1]):
                print(ctgn + ": " + str(contigcounts[ctgn]))

    def filter_reverse_small_contigs(self, size):
        for rid,read in self.lreads.items():
            toremove = []
            for item in read["maps"]:
                if item["ecc"] - item["scc"] < size and item["strand"] == 1:
                    toremove.append(item)
            for item in toremove:
                self.remove_contig_from_read(rid, item)

    
    def sort_by_starts(self):
        for read in self.lreads.values():
            read["maps"] = sorted(read["maps"], key = lambda x: x["scr"])

    def identify_hairpins(self):
        potential_hairpins = Counter()
        for rid, read in self.lreads.items():
            seen_contigs = {}
            for contig in read["maps"]:
                if contig["name"] in seen_contigs and seen_contigs[contig["name"]] != contig["strand"]:
                    potential_hairpins[rid] += 1
                if self.cellline in contig["name"]:
                    seen_contigs[contig["name"]] = contig["strand"]
        return potential_hairpins

    def get_problem_contigs(self):
        problem_contigs = defaultdict(set)
        for rid, read in self.lreads.items():
            seen_contigs = set()
            for contig in read["maps"]:
                if contig["name"] in seen_contigs:
                    problem_contigs[contig["name"]].add(rid)
                if self.cellline in contig["name"]:
                    seen_contigs.add(contig["name"])
        return problem_contigs

    # only use this when the maps have been sorted
    def get_maps_right(self, rid, start, end, cellline_only = True):
        maps = []
        for contig in self.lreads[rid]["maps"]:
            if contig["scr"] < start:
                continue
            elif contig["scr"] >= start and contig["scr"] <= end:
                if cellline_only and self.cellline not in contig["name"]:
                    continue
                else:
                    maps.append(contig)
            else: # contig["scr"] > end:
                return maps
        return maps # if less than 'end' bases to the right of the contig

    def get_overlapping_bases(self, ctg1, ctgs2, offset):
        overlapping = 0
        start = ctg1["scr"] + offset if ctg1["scr"] + offset > 0 else 0
        end = ctg1["ecr"] + offset
        #print(ctg1)
        #print("start: " + str(start))
        #print("end: " + str(end))
        for ctg2 in sorted(ctgs2, key= lambda x: x["scr"]):
            if ctg2["scr"] > end:
                return overlapping
            elif ctg2["ecr"] < start:
                continue
            else:
                if ctg2["ecr"] > end:
                    if ctg2["scr"] < start:
                        overlapping += end - start
                    elif ctg2["scr"] <= end:
                        overlapping += end - ctg2["scr"]
                else:
                    if ctg2["scr"] > start:
                        overlapping += ctg2["ecr"] - ctg2["scr"]
                    elif ctg2["ecr"] > start:
                        overlapping += ctg2["ecr"] - start
        return overlapping

    # pseudoalign matchin contigs if possible, punish overlapping mismatching contigs
    def pseudoalign(self,rid1, rid2, offset, debug= False):
        #print("-"*50)
        #print("offset: " + str(offset))
        overall_score = 0
        for ctg1 in sorted(self.lreads[rid1]["maps"], key=lambda x: x["scr"]):
            #print("\t".join([ctg1["name"], str(ctg1["scr"]+offset), str(ctg1["ecr"]+offset)]))
            #print(self.lreads[rid2]["maps"][0])
            #print(self.lreads[rid2]["maps"][-1])
            if ctg1["ecr"] + offset < self.lreads[rid2]["maps"][0]["scr"] or ctg1["scr"] + offset > self.lreads[rid2]["maps"][-1]["ecr"]:
               
               continue
            best_pair_score = -1
            for ctg2 in self.lreads[rid2]["mapsc"][ctg1["name"]]:
                vstart1 = ctg1["scr"] - ctg1["scc"] 
                vstart2 = ctg2["scr"] - ctg2["scc"] 
                distance = vstart2 - vstart1 - offset
                scoret = -abs(distance)
                scoret += 2*min(ctg1["ecc"]-ctg1["scc"], ctg2["ecc"] - ctg2["scc"])
                #print("\t".join([ctg1["name"], ctg2["name"], "scoret", str(scoret),"distance",str(distance)]))
                if scoret > best_pair_score:
                    best_pair_score = scoret 
            #print("\t".join(["best_pair_score", ctg1["name"], str(best_pair_score)]))
            if best_pair_score >= 0:
                overall_score += best_pair_score
            else: # check if problem for that contig
                #print(ctg1["name"])
                #print("getting overlapping bases for: " + str(ctg1["name"]))
                nscore = - 3*self.get_overlapping_bases(ctg1, self.lreads[rid2]["maps"], offset)
                #print("penalty: " + str(nscore))
                overall_score += nscore
        #print("overall_score: " + str(overall_score))
        return overall_score

    def get_possible_offsets(self, lr1, lr2, min1=0):
        """ Gets all possible offsets for two specific reads given its contig signatures
        min1 can be set to exclude all contigs with coordinates lower than that value
        """
        samectgs = self.lreads[lr1]["ctgset"] & self.lreads[lr2]["ctgset"]
        poffs = []
        if not samectgs:
            return []
        else:
            for ctgn in samectgs:
                ctgs1 = self.lreads[lr1]["mapsc"][ctgn]
                ctgs2 = self.lreads[lr2]["mapsc"][ctgn]
                for ctg1, ctg2 in product(ctgs1, ctgs2):
                    d1 = ctg1["scr"] - ctg1["scc"]
                    if d1 < min1:
                        continue
                    d2 = ctg2["scr"] - ctg2["scc"]
                    poffs.append((d1,d2))
            roffs = []
            for off in poffs:
                for coff in roffs:
                    if abs( (off[1]-off[0]) - (coff[1]-coff[0]) ) < 100:
                        break
                else:
                    roffs.append(off)
            return roffs


    # Returns the space in y that it needs, depends on the number of longreads that this object contains
    def to_SVG(self, svglongread, lread_ids, ctglengths, xoff, yoff, ctg_y_drawsize=12, show_lr_ids=False):
        zoom = svglongread.zoom
        img = svglongread.dwg
        ypad = 7
        xpad = 20
        col = "black"
        nr_longreads = len(self.lreads)
        ypos = yoff
        # TODO reimplement
        #if show_lr_ids:
            #for lrid, lrc in self.longread_coords.items():
                #ylen = nr_longreads * y_space_per_lr - ypos
                #rect = img.add(svgwrite.shapes.Rect((xoff+(lrc[0]/100),yoff+ypos), ((lrc[1]-lrc[0])/100,ylen+ypad), stroke='green', stroke_width=1 ))
                #rect.fill(color="none").dasharray([2, 2])
                #img.add(img.text(lrid, insert=(xoff+(lrc[0]/100),yoff+ypos-1),fill="green", style="font-size:2"))
                #ypos += y_space_per_lr
        ypos += ctg_y_drawsize + ypad
        ctg_y_halfdrawsize = ctg_y_drawsize/2
        ctg_relative_positions = cycle([-ctg_y_halfdrawsize-2, ctg_y_halfdrawsize+4, -ctg_y_halfdrawsize-5, ctg_y_halfdrawsize+7])
        gradient_idc = 0

        def get_colors_from_nr(nr):
            col1 = "#000000"
            col2 = "#000000"
            if nr==0:
                col1 = "#FF0000"
                col2 = "#0000FF"
            elif nr == 1:
                col1 = "#5CDB95"
                col2 = "#05386B"
            elif nr==2:
                col1 = "#FFE119"
                col2 = "#000000"
            elif nr == 3:
                col1 = "#379683"
                col2 = "#BC986A"
            elif nr==4:
                col1 = "#911eb4"
                col2 = "#a9a9a9"
            elif nr == 5:
                col1 = "#659DBD"
                col2 = "#F64C72"
            elif nr==6:
                col1 = "#000075"
                col2 = "#f58231"
            elif nr == 7:
                col1 = "#553D67"
                col2 = "#C38D9E"
            elif nr==8 :
                col1 = "#e6beff"
                col2 = "#808000"
            elif nr == 9:
                col1 = "#41B3A3"
                col2 = "#501B1D"
            return [col1, col2]


        def get_colors(ctgn):
            col1, col2 = get_colors_from_nr(int(ctgn.rstrip(self.cellline)[-1]))
            total = 0
            return [col1, col2]

        for lrid in lread_ids:
            lread = self.lreads[lrid]
            img.add(img.text(lrid, insert = (xoff - 50, ypos), style="font-size:5"))
            g = img.defs.add(Group(id=lrid))
            g.add(svgwrite.shapes.Line((xoff, ypos), ( xoff + lread["length"]/zoom, ypos), stroke=svgwrite.rgb(0, 0, 0, '%')))
            for contig in sorted(lread["maps"], key= lambda x: x["scr"]):
                #print(read)
                sc = contig["scr"]
                ec = contig["ecr"]
                scc = contig["scc"]
                ecc = contig["ecc"]
                ctgn = contig["name"]
                #ctg = read[0]
                if ctgn.startswith("chr"):
                    ctgnd = ctgn[0:9]
                else:
                    ctgnd = "$" + ctgn + "q"

                gradient_idc += 1
                gradient_id1 = str(gradient_idc)
                lineargrad1 = img.defs.add(svgwrite.gradients.LinearGradient(id=gradient_id1 , x1=-scc/(ecc-scc), x2=1+(ctglengths[shortname(ctgn)]-ecc)/(ecc-scc), y1=0, y2=0))
                col1, col2 = get_colors(shortname(ctgn))
                
                lineargrad1.add_stop_color("0%",col1)
                lineargrad1.add_stop_color("50%","#FFFFFF")
                lineargrad1.add_stop_color("100%",col2)

                g.add(svgwrite.shapes.Rect((xoff+sc/zoom,ypos-ctg_y_halfdrawsize), ((ec-sc)/zoom,ctg_y_drawsize), stroke='black', stroke_width=1, fill='url(#'+str(gradient_idc)+')'))
                yt = ypos + next(ctg_relative_positions)
                g.add(img.text(ctgnd, insert=(xoff+(sc/zoom),yt),fill=col, style="font-size:3"))
                if contig["strand"] == 0:
                    direction = ">"
                else:
                    direction = "<"
                g.add(img.text(direction, insert=(xoff+sc/zoom,ypos+2),style="font-size:6"))
            img.add(g)
            ypos += ctg_y_drawsize + 30
        return ypos-yoff

    
    def pseudoalign_all(self, debug=False):
        dists = defaultdict(dict)
        lr_scores = defaultdict(dict)
        toalign = OrderedDict(sorted(self.lreads.items(), key= lambda lr: lr[1]["length"]))
        while toalign:
            lr1, lread1 = toalign.popitem()
            #print(lread1["length"])
            dists[lr1][lr1] = 0
            lr_scores[lr1][lr1] = 0
            for lr2 in toalign:
                if lr2 == lr1:
                    continue
                offs = self.get_possible_offsets(lr1,lr2)
                scores = []
                for d1, d2 in offs:
                    scores.append(self.pseudoalign(lr1,lr2, d2-d1))
                if scores:
                    if debug:
                        print("-"*40)
                        print(lr1 + "\t" + lr2)
                        print(offs)
                        print(scores)
                    if max(scores) > 0:
                        sidx = scores.index(max(scores))
                        dists[lr1][lr2] = offs[sidx][0]-offs[sidx][1] # save offset in table, score doesn't matter 
                        dists[lr2][lr1] = offs[sidx][1]-offs[sidx][0] # save offset in table, score doesn't matter 
                    elif max(scores) == 0:
                        dists[lr1][lr2] = None # distance can't be determined, but no contradiction
                        dists[lr2][lr1] = None # distance can't be determined, but no contradiction
                    lr_scores[lr1][lr2] = max(scores)
                    lr_scores[lr2][lr1] = max(scores)
        return (lr_scores, dists)




    def remove_problem_contigs(self,ctgn, prids):
        #print(self.ctg2lreads[ctg])
        print("contig " + ctgn)
        print("prids " + str(prids))

        #prids has all reads which have the contig doubly
        #let's see if evidence can be found, that both make sense
        singles = self.ctg2lreads[ctgn] - set(prids)
        print("singles: " + str(singles))
        for prid in prids:
            ctgpos1 = []
            ctgs1 = []
            for ctg in self.lreads[prid]["mapsc"][ctgn]:
                ctgpos1.append(ctg["scr"] - ctg["scc"])
                ctgs1.append(ctg)
            vals = {}
            for d1 in ctgpos1:
                vals[d1] = []
                for rid in singles:
                    ctg2 = self.lreads[rid]["mapsc"][ctgn][0]
                    d2 = ctg2["scr"] - ctg2["scc"]
                    offset = d2 - d1
                    val = self.pseudoalign(prid, rid, offset)
                    vals[d1].append(val)
            for idx, d1 in enumerate(ctgpos1):
                if max(vals[d1]) < 0:
                    print("removing contig " + ctgs1[idx]["name"] + " from " + prid)
                    self.remove_contig_from_read(prid, ctgs1[idx])
                else:
                    print(vals)
                    
            sys.exit()

            #print(vals)


    def sort_contigs_in_reads(self):
        # sort contigs by left coordinate
        for rid in self.lreads:
            self.lreads[rid]["maps"] = sorted(self.lreads[rid]["maps"], key = lambda x: x["scr"])

    def turn_longreads_around(self, revs=[], logging=False):
        for rid in self.lreads:
            fcount = 0
            for contig in self.lreads[rid]["maps"]:
                if contig["name"] in revs:
                    fcount += (contig["strand"]*2 - 1) # +1 if forward, -1 if revcomp
                else:
                    fcount += (contig["strand"]*(-2) + 1) # -1 if forward, +1 if revcomp
            if fcount < 0:
                if logging:
                    logging.write("Turning around: " + rid + "\n")
                length = self.lreads[rid]["length"]
                self.lreads[rid]["reverse"] = True
                for contig in self.lreads[rid]["maps"]:
                    tmp = contig["scr"]
                    contig["scr"] = length - contig["ecr"] 
                    contig["ecr"] = length - tmp
                    contig["strand"] = 0 if contig["strand"] == 1 else 1


class LongReadSVG(object):

    def __init__(self, filehandle, zoom=100,shortread_info=False):
        self.zoom = zoom
        dwg = svgwrite.Drawing(filehandle, size=(u'1700', u'2200'), profile='full')
           
        #svgwrite.gradients.LinearGradient((0,0), (100,0))
        yp = 10
        self.xtext = 10
        xpad = 20
        self.xpad = xpad
        dwg.add(dwg.text("600 bases", insert=(self.xpad, yp), fill='black', style="font-size:7"))
        yp += 5 
        dwg.add(dwg.line((self.xpad, yp), ( xpad + 600/zoom, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
        dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
        dwg.add(dwg.line((xpad + 600/zoom, yp-2), ( xpad +600/zoom, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
        yp += 10
        dwg.add(dwg.text("10000 bases", insert=( xpad, yp), fill='black', style="font-size:7"))
        yp += 5
        dwg.add(dwg.line((xpad, yp), ( xpad + 10000/zoom, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
        dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
        for i in range(1,10):
            dwg.add(dwg.line( (xpad + (10000/zoom)/10 * i, yp-1), (xpad + (10000/zoom)/10 * i, yp+1), stroke=svgwrite.rgb(0,0,0,'%')))
        dwg.add(dwg.line((xpad + 10000/zoom, yp-2), ( xpad +10000/zoom, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
        yp += 10
        dwg.add(dwg.text("100k bases", insert=( xpad, yp), fill='black', style="font-size:7"))
        yp += 5
        dwg.add(dwg.line((xpad, yp), ( xpad + 100000/zoom, yp), stroke=svgwrite.rgb(0, 0, 0, '%')))
        dwg.add(dwg.line((xpad, yp-2), ( xpad , yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
        for i in range(1,10):
            dwg.add(dwg.line( (xpad + (100000/zoom)/10 * i, yp-1), (xpad + (100000/zoom)/10 * i, yp+1), stroke=svgwrite.rgb(0,0,0,'%')))
        dwg.add(dwg.line((xpad + 100000/zoom, yp-2), ( xpad +100000/zoom, yp+2), stroke=svgwrite.rgb(0, 0, 0, '%')))
        yp += 10
        g1M = dwg.defs.add(dwg.g(id='g003'))
        g1M.add(dwg.text("M bases", insert=( xpad, yp), fill='black', style="font-size:7"))
        yp += 5
        g1M.add(dwg.line((xpad, yp), ( xpad + 5000000/zoom, yp), stroke=svgwrite.rgb(0, 0, 0, '%'), stroke_width="5"))
        g1M.add(dwg.line((xpad, yp-7), ( xpad , yp+7), stroke=svgwrite.rgb(0, 0, 0, '%'), stroke_width="5"))
        for i in range(0,51):
            if i%10 == 0 and i != 0:
                g1M.add(dwg.line( (xpad + (1000000/zoom)/10 * i, yp-7), (xpad + (1000000/zoom)/10 * i, yp+200), stroke=svgwrite.rgb(0,0,0,'%'),stroke_width="7"))
            else:
                g1M.add(dwg.line( (xpad + (1000000/zoom)/10 * i, yp-5), (xpad + (1000000/zoom)/10 * i, yp+200), stroke=svgwrite.rgb(0,0,0,'%'),stroke_width="2"))
        g1M.add(dwg.text("1 M", insert=( xpad+ 1000000/zoom-10, yp-10), fill='black', style="font-size:100"))
        g1M.add(dwg.text("2 M", insert=( xpad+ 2000000/zoom-10, yp-10), fill='black', style="font-size:100"))
        g1M.add(dwg.text("3 M", insert=( xpad+ 3000000/zoom-10, yp-10), fill='black', style="font-size:100"))
        g1M.add(dwg.text("4 M", insert=( xpad+ 4000000/zoom-10, yp-10), fill='black', style="font-size:100"))
        g1M.add(dwg.text("5 M", insert=( xpad+ 5000000/zoom-10, yp-10), fill='black', style="font-size:100"))
        dwg.add(g1M)
        #yp += 20
        #rect = dwg.add(svgwrite.shapes.Rect((xpad,yp-3), (2000/zoom,7), stroke='green', stroke_width=1 ))
        #rect.fill(color="none").dasharray([2, 2])
        #dwg.add(dwg.text("region of scaffold covered by long read", insert=( xpad+ 2000/zoom+ 5, yp+2), fill='black', style="font-size:7"))
        yp += 10
        if shortread_info:
            dwg.add(svgwrite.shapes.Rect((xpad,yp-3), (2000/zoom,7), stroke='gray', stroke_width=1, fill='none' ))
            dwg.add(dwg.text("contigs from short read data", insert=( xpad+ 2000/zoom + 5, yp+2), fill='black', style="font-size:7"))
            yp += 10
        dwg.add(svgwrite.shapes.Rect((xpad,yp-3), (2000/zoom,7), stroke='black', stroke_width=1, fill='none' ))
        dwg.add(dwg.text("contigs from long read data", insert=( xpad+ 2000/zoom + 5, yp+2), fill='black', style="font-size:7"))
        yp += 20
        self.yp = yp
        self.dwg = dwg
