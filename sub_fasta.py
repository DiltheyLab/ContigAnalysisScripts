from argparse import ArgumentParser
from Bio import SeqIO
from scaffold import revcomp



parser = ArgumentParser()
parser.add_argument("inputfile", help="Input File in FASTA format")
parser.add_argument("chromosom", help="Identifier")
parser.add_argument("start", type=int, help="Start coordinate")
parser.add_argument("end", type=int, help="End coordinate")
args = parser.parse_args()

chroms = {}

for read in SeqIO.parse(args.inputfile, "fasta"):
    chroms[read.id] = str(read.seq)

for chrom in chroms:
    print(chrom)

print(chroms[args.chromosom][args.start:args.end])
print(revcomp(chroms[args.chromosom])[args.start:args.end])
