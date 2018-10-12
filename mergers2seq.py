from argparse import ArgumentParser
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
import numpy as np
import graphviz as gv
import networkx as nx


parser = ArgumentParser()
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

with open(args.input) as f:
    for line in f:
        [modestring, node1, l1, node2, l2, distance, newnode] = line.split()
        if not node1.startswith("cluster"):
            dot.add_node(node1, length = l1)
        if not node2.startswith("cluster"):
            dot.add_node(node2, length = l2)
        dot.add_node(newnode, mode = modestring)
        dot.add_edge(node1, newnode)
        dot.add_edge(newnode, node1, link = "left")
        dot.add_edge(node2,newnode)
        if modestring == "extension":
            new_length = distance + l2
            dot.add_edge(newnode, node2, link = "right")
        else:
            new_length = l1
        dot.nodes[newnode]["length"] = new_length
        

#print(dot.nodes["cluster_1"])
#print(dot.nodes["653b2764-e21e-458a-99e4-702f7eac052a"])

# delete a node
"""
node = "cluster_101"
parents = list(dot.predecessors(node))
children = list(dot.successors(node))
child = children[0]
print(parents)
print(children)
#for parent in parents:
#    dot[parent]
dot.remove_node(node)
for parent in parents:
    print(parent)
    dot.add_edge(parent, child)
    for item in dot.successors(parent):
        print("child: " + str(item))
"""



def get_start_node(g, end_node):
    n = end_node
    while True:
        print("getting start node from: " + n)
        for node2, edge in g[n].items():
            if "link" in edge and edge["link"] == "left":
                n = node2
                break
        else:
            return n


# get final nodes
final_nodes = []
for node in dot.nodes:
    for edge in dot[node]:
        if "link" in dot[node][edge]:
            continue
        break
    else:
        final_nodes.append(node)

    
def suffix(path,offset):
    suffix = []
    for part, name in path:
        if part[0] <= offset and part[1] > offset:
            new_path.append((part[0],offset), name,(ctgs,ctge))
            return new_path
        elif part[0] > offset:
            suffix.append(((part[0]-offset, part[1]-offset), name))
    return suffix

# move to next extension node
# up until n2 !!
def move_to_extension(g, n1, n2):
    print("Move to Extension " + n1 + " and " + n2)
    current_node = go_down(g, n1)
    #print("current_node: " + current_node)

    while True:
        #print(current_node + " is  current_node" )
        if current_node == n2:
            return current_node
        if "mode" in g.nodes[current_node] and g.nodes[current_node]["mode"] == "extension":
            print("extension found at node: " + current_node)
            return current_node
        else:
            current_node = go_down(g, current_node)
            #print("current_node: " + current_node)

def go_right(g, n):
    for n2 in g[n]:
        if "link" in g[n][n2] and g[n][n2]["link"] == "right":
            return n2
    else:
        print("Problem finding right neigbour")
        sys.exit(0)

def go_down(g, n):
    #print("going down from: " + n)
    for n2 in g[n]:
        if "link" in g[n][n2]:
            continue
        else:
            return n2
    else:
        return n2
                
                
def find_path(g, n1, n2):
    print("Finding path between " + n1 + " and " + n2)
    path = [((1,g.nodes[n1]["length"]),n1, 0)]
    if n1 == n2:
        return path
    current_node = n1
    while True:
        current_node = move_to_extension(g, current_node, n2)
        if "mode" in g.nodes[current_node] and g.nodes[current_node]["mode"] == "extension":
            rn = go_right(g, current_node)
            path.append(find_path(g, get_start_node(g, rn) ,rn))
        if current_node == n2:
            return path

#print(dot.nodes)
print(dot.nodes["cluster_1"]["mode"])
print(final_nodes[0])
sn = get_start_node(dot, final_nodes[0])
print(dot["cluster_277"]["cluster_276"])
print(dot.nodes["cluster_277"])
print(dot.nodes["cluster_276"])
#sys.exit(0)
print(find_path(dot, sn, final_nodes[0]))

# TODO
'''
find_path n1 n2:
    while not n2:
    move to extension node e from n1
    extension = suffix( find_path(get_start_node(e),e))
    return n1
'''




# translate to dot and plot with graphviz
dotgv = gv.Digraph(comment="clusters")
for node in dot:
    #print(node)
    #print(dot[node])
    if "mode" in (dot.nodes[node]) and dot.nodes[node]["mode"] == "extension":
        dotgv.node(node, color = "red", style = "filled")
    else:
        dotgv.node(node)

    for edge in dot[node]:
        dotgv.edge(node,edge)
    
dotgv.format = "svg"
dotgv.render('clusters_reduced.gv', view=False)  # doctest: +SKIP




#for item in dot:
#    for edge in dot[item]:
#        print(dot[item][edge])

#print(nx.DiGraph.predecessors(dot["cluster_101"]))


#for item in dot:
#    print(item)



