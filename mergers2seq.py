from argparse import ArgumentParser
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
import numpy as np
from graphviz import Digraph


parser = ArgumentParser()
parser.add_argument("input", help="fastq-inputfile")
parser.add_argument("outfile", help="fasta-outputfile for scaffolds")

args = parser.parse_args()

""" thaz how it look like
incorporation	cluster_6	d625a932-ce76-4d38-a486-c708414e4c37	193950	cluster_7
incorporation	cluster_7	ebb2adf5-7217-4f60-aca3-1f3e6205ee41	211867	cluster_8
incorporation	cluster_8	02fa0a3d-0721-4017-948f-3eeff063c6c9	203035	cluster_9
extension	cluster_9	ab0cb3e0-e72a-4ac9-8d3d-f86cc2867a7e	239134	cluster_10
"""


dot = Digraph(comment='Them clusters')

with open(args.input) as f:
    for line in f:
        [mode, node1, node2, distance, newnode] = line.split()
        if not node1.startswith("cluster"):
            dot.node(node1, node1)
        if not node2.startswith("cluster"):
            dot.node(node2, node2)
        dot.node(newnode, newnode)
        dot.edge(node1,newnode)
        dot.edge(node2,newnode)


print(dot.source)  # doctest: +NORMALIZE_WHITESPACE

dot.render('clusters.gv', view=True)  # doctest: +SKIP



