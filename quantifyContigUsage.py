from Bio import SeqIO

contigs = {}

for read in SeqIO.parse("QBLContigs.fasta", "fasta"):
    contigs[read.id] = len(read.seq)

tsum = 0
for contig in contigs:
    tsum += contigs[contig]

kcontigs = set()


with open("all.mmi.0x800.erates") as f:
    for line in f:
        sline = line.split()
        ctg = sline[1]
        if ctg.endswith("QBL"):
            if not ctg in kcontigs:
                kcontigs.add(ctg)
        

ssum = 0
for contig in kcontigs:
    ssum += contigs[contig]

print("total bases in contigs: " + str(tsum))
print("bases in contigs: " + str(ssum))
    
