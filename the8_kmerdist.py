from argparse import ArgumentParser
from Bio import SeqIO
from bitstring import BitArray

parser = ArgumentParser()
parser.add_argument("--k", type= int, default = 100)

args = parser.parse_args()
k = args.k

celllines = {"MOU": BitArray(bin='00000001'), "QBL": BitArray(bin='00000010'), "PGF": BitArray(bin='00000100'), "APD": BitArray(bin='00001000'), "SSTO": BitArray(bin='00010000'), "MCF": BitArray(bin='00100000'), "DBB": BitArray(bin='01000000')}
#celllines = {"MOU_test": 0b00000001, "QBL_test": 0b00000010}
#celllines = {"MOU_test": BitArray(bin='00000001'), "QBL_test": BitArray(bin='00000010')}
#celllines = {"MOU_test": BitArray(bin='00000001'), "QBL_test": BitArray(bin='00000010')}

counts = {1: 0, 2:0,3:0,4:0,5:0,6:0,7:0}

kmers = {}

def set_kmer(cellline, kmer):
    if not kmer in kmers:
        kmers[kmer] = BitArray(bin='00000000')
    kmers[kmer] |= celllines[cellline]

for cline in celllines.keys():
    print("Analyzing " + cline + " ...")
    filename = cline + "Contigs.fasta"
    for read in SeqIO.parse(filename, "fasta"):
        for i in range(k, len(read.seq)):
            kmer = str(read.seq)[i-k:i]
            kmerrc = str(read.reverse_complement().seq)[i-k:i]
            set_kmer(cline, kmer)
            set_kmer(cline, kmerrc)


for kmer in kmers:
    counts[kmers[kmer].count(1)] += 1

print(counts)
