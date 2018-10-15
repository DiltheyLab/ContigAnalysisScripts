from argparse import ArgumentParser
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
import numpy as np
import graphviz as gv
import networkx as nx


parser = ArgumentParser()
parser.add_argument("mergefile", help="merge file, output of less_naive_scaffolding.py")
parser.add_argument("input", help="fastq-inputfile")
parser.add_argument("outfile", help="fasta-outputfile for scaffolds")

args = parser.parse_args()

""" thaz how it look like
incorporation	cluster_6	280489	d625a932-ce76-4d38-a486-c708414e4c37	19987	193950	cluster_7
incorporation	cluster_7	280489	ebb2adf5-7217-4f60-aca3-1f3e6205ee41	3214	211867	cluster_8
incorporation	cluster_8	280489	02fa0a3d-0721-4017-948f-3eeff063c6c9	6013	203035	cluster_9
extension	cluster_9	280489	ab0cb3e0-e72a-4ac9-8d3d-f86cc2867a7e	104531	239134	cluster_10

"""


dot = nx.DiGraph()

with open(args.mergefile) as f:
    for line in f:
        [modestring, node1, l1, node2, l2, distance, newnode] = line.split()
        if not node1.startswith("cluster"):
            dot.add_node(node1)
            dot.nodes[node1]["length"] = int(l1)
        if not node2.startswith("cluster"):
            dot.add_node(node2)
            dot.nodes[node2]["length"] = int(l2)
        dot.add_node(newnode, mode = modestring)
        dot.add_edge(node1, newnode, link = "down")
        dot.add_edge(newnode, node1, link = "left")
        dot.add_edge(node2, newnode, link = "down")
        dot.add_edge(newnode, node2, link = "right")
        if modestring == "extension":
            new_length = int(distance) + int(l2)
        else:
            new_length = int(l1)
        dot.nodes[newnode]["length"] = new_length
        
print("cluster_660 len: " + str(dot.nodes["cluster_660"]["length"]))
print("8a3dd len: " + str(dot.nodes["8a3dd374-8c36-4ceb-88ac-ff02534856d0"]["length"]))

def get_start_node(g, end_node):
    n = end_node
    while True:
        #print("getting start node from: " + n)
        for node2, edge in g[n].items():
            if edge["link"] == "left":
                n = node2
                break
        else:
            return n

# move to next extension node
# up until n2 !!
def move_to_extension(g, n1, n2):
    current_node = n1
    while current_node != n2:
        current_node = go_down(g, current_node)
        #print(current_node + " is  current_node" )
        if "mode" in g.nodes[current_node] and g.nodes[current_node]["mode"] == "extension":
            #print("extension found at node: " + current_node)
            return current_node
            #print("current_node: " + current_node)
    return False

def go_right(g, n):
    for n2 in g[n]:
        if g[n][n2]["link"] == "right":
            return n2
    else:
        print("Problem finding right neigbour")
        sys.exit(0)

def go_left(g, n):
    for n2 in g[n]:
        if g[n][n2]["link"] == "left":
            return n2
    else:
        print("Problem finding left neigbour")
        sys.exit(0)

def go_down(g, n):
    for n2 in g[n]:
        if g[n][n2]["link"] == "down":
            return n2
    return n

def get_suffix(path, dist, offset):
    npath = []
    print("path to be suffixed: " + str(path))
    print("offset: " + str(offset))
    print("dist: " + str(dist))
    for segment in path:
        coords, name, seg_offset = segment
        if int(coords[0]) + int(dist) > int(offset):
            npath.append(((int(coords[0]) + int(dist) + 1, int(coords[1]) + int(dist) + 1), name, int(seg_offset)))
        elif int(coords[0]) + int(dist) <= int(offset) and int(coords[1]) + int(dist) >= int(offset):
            npath.append(((int(offset) + 1, int(dist) + int(coords[1]) ), name, int(seg_offset) + int(dist)))
    return npath

                
def find_path(g, n1, n2):
    #print("Finding path between " + n1 + " and " + n2)
    path = [((0, int(g.nodes[n1]["length"] )), n1, 0)]
    if n1 == n2:
        return path
    current_node = n1
    while True:
        a =  move_to_extension(g, current_node, n2)
        if a:
            current_node = a
            rightn = go_right(g, current_node)
            leftn = go_left(g, current_node)
            nl = g.nodes[current_node]["length"] 
            ol = g.nodes[leftn]["length"]
            distance = int(nl) - int(ol)
            assert distance > 0
            # here one could easily introduce more sophisticated merging via pairwise alignment
            extension_path = get_suffix(find_path(g, get_start_node(g, rightn), rightn), distance, ol)
            path.extend(extension_path)
        else:
            return path


# get final nodes
final_nodes = []
for node in dot.nodes:
    #print("Getting final node for: " + node)
    n1 = node
    n2 = go_down(dot, n1)
    while n1 != n2:
        n1 = n2
        n2 = go_down(dot, n2)
    if n1 not in final_nodes:
        final_nodes.append(n1)

for node in final_nodes:
    sn = get_start_node(dot, node)
    dot.nodes[sn]["color"] = "green"

#print(final_nodes)

easy_ones = ["cluster_544", "cluster_543", "cluster_701"]
hard = ["cluster_697"]
hard = ["cluster_706"]
easy = ["cluster_660"]
paths = {}
for i in easy:
    #fn = "cluster_" + str(i)
    fn = i
    print("final node: " + fn)
    sn = get_start_node(dot, fn)
    print("starts at: " + sn)
    #sys.exit(0)
    solpath = find_path(dot, sn, fn)
    print(solpath)
    paths[fn] = solpath

#sys.exit(0)

print("Loading sequences ... ") 
lrs = {}
with open(args.input) as f:
    for i, line in enumerate(f):
        if i % 4 == 0:
            lrid = line.split(" ")[0][1:]
        elif i % 4 == 1:
            lrs[lrid] = line.rstrip()
    
print("done")



# on to the meat
with open(args.outfile, "w") as f:
    for i, path in paths.items():
        f.write(">" + i + "\n")
        for segment in path:
            f.write(lrs[segment[1]][segment[2]:])
        f.write("\n")

# translate to dot and plot with graphviz
dotgv = gv.Digraph(comment="clusters")
for node in dot:
    #print(node)
    #print(dot[node])
    if "mode" in (dot.nodes[node]) and dot.nodes[node]["mode"] == "extension":
        dotgv.node(node, color = "red", style = "filled")
    elif "color" in dot.nodes[node]:
        dotgv.node(node, color = dot.nodes[node]["color"], style = "filled")
    else:
        dotgv.node(node)

    for edge in dot[node]:
        if dot[edge][node]["link"] == "left":
            dotgv.edge(node,edge, "l")
        elif dot[node][edge]["link"] == "down":
            dotgv.edge(node,edge)
    
dotgv.format = "svg"
dotgv.render('clusters_reduced.gv', view=False)  # doctest: +SKIP


#for item in dot:
#    for edge in dot[item]:
#        print(dot[item][edge])

#print(nx.DiGraph.predecessors(dot["cluster_101"]))


#for item in dot:
#    print(item)



