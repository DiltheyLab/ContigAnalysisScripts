from Bio import SeqIO
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("ContigFile")
parser.add_argument("ErateFile")
parser.add_argument("linename")
args = parser.parse_args()

contigs = {}
for read in SeqIO.parse(args.ContigFile, "fasta"):
    contigs[read.id] = len(read.seq)

tsum = 0
for contig in contigs:
    tsum += contigs[contig]

kcontigs = set()


with open(args.ErateFile) as f:
    for line in f:
        sline = line.split()
        ctg = sline[1]
        if ctg.endswith(args.linename):
            if not ctg in kcontigs:
                kcontigs.add(ctg)
        

ssum = 0
for contig in kcontigs:
    ssum += contigs[contig]

print("total bases in contigs: " + str(tsum))
print("bases in contigs: " + str(ssum))
print("percentage: " + str(ssum/tsum))
    
