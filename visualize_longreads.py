from argparse import ArgumentParser
import pickle
import sys
import matplotlib.pyplot as plt
import matplotlib
import random 
from collections import defaultdict, deque
import squarify
from scaffold import Scaffold, Longreads, LongReadSVG

parser = ArgumentParser()
parser.add_argument("inputfile", help="Input Pickle File with LongReads.")
parser.add_argument("squareplot", help="Plots Long Read lengths in squarify plot")
args = parser.parse_args()

with open(args.inputfile, "rb") as f:
    scafs = pickle.load(f)

def stringify(nr):
    if nr >= 1000000:
        return "{:.1f}".format(nr/1000000) + " M"
    elif nr >= 10000:
        return "{:.0f}".format(nr/1000) + " k"
    elif nr > 1000:
        return "{:.1f}".format(nr/1000) + " k"
    else:
        return str(nr)

lr_lengths = []
for lr in scafs.lreads.values():
    lr_lengths.append(lr["length"])
norm = matplotlib.colors.Normalize(vmin=min(lr_lengths), vmax=max(lr_lengths))
lr_colors = [(random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for x in range(len(lr_lengths))]
#rs = [(random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for x in range(len(lr_lengths))
#matplotlib.cm.gist_ncar(value) for value in random.sample(norm,len(lr_lengths))]


plt.rc('font', size=15)          # controls default text sizes
#plt.subplot(121)
squarify.plot(sizes=lr_lengths, label=[stringify(x) for x in lr_lengths], alpha=.9,color = lr_colors )
plt.axis('off')
plt.savefig(args.squareplot)
