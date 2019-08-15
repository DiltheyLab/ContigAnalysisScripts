from argparse import ArgumentParser
from Bio import SeqIO
import sys


parser = ArgumentParser()
parser.add_argument("inputfile", help="Input Files in FASTA format")
parser.add_argument("id", help="ID of read in FASTA file")
parser.add_argument("--start", type=int, help="Starting coordinate")
parser.add_argument("--end", type=int, help="End coordinate")
parser.add_argument("output", help="Output File in FASTA format")

args = parser.parse_args()

sequence = {}
for read in SeqIO.parse(args.inputfile, "fasta"):
    if args.id == read.id:
        sequence[args.id] = str(read.seq)
        break
else:
    print("id not found")
    sys.exit()

start = args.start if args.start else 0
end = args.end+1 if args.end else len(sequence[args.id])+1

with open(args.output, "w") as f:
    nid = ">" + args.id + "_subseq_" + str(start) + ":" + str(end)
    f.write(nid + "\n")
    f.write(sequence[args.id][start:end])


