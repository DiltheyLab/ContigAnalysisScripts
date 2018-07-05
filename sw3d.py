import argparse
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from bidict import bidict


parser = argparse.ArgumentParser(description = "Calculates the multiple sequence alignment for three sequences using dynamic programming. ONLY calculates the MSA for the first sequence in each input file.") 
parser.add_argument("input1", help = "Shold be the long read (Nanopore, Pacbio)")
parser.add_argument("input2", help = "Contig 1")
parser.add_argument("input3", help = "Contig 2")

args=parser.parse_args()

records = { }
rec = SeqIO.parse(args.input1, "fasta").__next__()
records[rec.id] = str(rec.seq.upper())
rec1 = rec.id
rec = SeqIO.parse(args.input2, "fasta").__next__()
records[rec.id] = str(rec.seq.upper())
rec2 = rec.id
rec = SeqIO.parse(args.input3, "fasta").__next__()
records[rec.id] = str(rec.seq.upper())
rec3 = rec.id

mat = np.zeros((len(records[rec1])+1,len(records[rec2])+1,len(records[rec3])+1), dtype=np.uint8)
traceback = np.zeros((len(records[rec1])+1,len(records[rec2])+1,len(records[rec3])+1), dtype=np.uint8)
gap = -1
match = 3 
lrmissmatch = -1
srmissmatch = -3

bits = bidict({'sx':0, 'sy':128, 'sz':256, 'x':1, 'y':2, 'z':4, 'xy': 8, 'xz': 16, 'yz': 32, 'xyz': 64})
moves = {'sx':(0,0,0), 'sy':(0,0,0), 'sz':(0,0,0),  'x':(1,0,0), 'y':(0,1,0), 'z':(0,0,1), 'xy': (1,1,0), 'xz': (1,0,1), 'yz': (0,1,1), 'xyz': (1,1,1)}
print("Calculating matrix...")

# Initialization of the three sides
for i in range(1, len(records[rec1])+1):
    l1 = records[rec1][i-1]
    for j in range(1, len(records[rec2])+1):
        l2 = records[rec2][j-1]
        s12 = match if l1 == l2 else lrmissmatch
        d = {'s':0, 'xy': mat[i-1][j-1][0] + s12, 'x': mat[i-1][j][0] + gap, 'y': mat[i][j-1][0] + gap}
        a = max(d, key=d.get)
        mat[i][j][0] = d[a]
for i in range(1, len(records[rec1])+1):
    l1 = records[rec1][i-1]
    for k in range(1, len(records[rec3])+1):
        l3 = records[rec3][k-1]
        s13 = match if l1 == l3 else lrmissmatch
        d = {'s':0, 'xz': mat[i-1][0][k-1] + s13, 'x': mat[i-1][0][k] + gap, 'z': mat[i][0][k-1] + gap}
        a = max(d, key=d.get)
        mat[i][0][k] = d[a]
for j in range(1, len(records[rec2])+1):
    l2 = records[rec2][j-1]
    for k in range(1, len(records[rec3])+1):
        l3 = records[rec3][k-1]
        s23 = match if l2 == l3 else srmissmatch
        d = {'s':0, 'yz': mat[0][j-1][k-1] + s23, 'y': mat[0][j-1][0] + gap, 'z': mat[0][j][k-1] + gap}
        a = max(d, key=d.get)
        mat[0][j][k] = d[a]



for i in range(1,len(records[rec1])+1):
    print("step " + str(i) + "/" + str(len(records[rec1])))
    l1 = records[rec1][i-1]
    for j in range(1,len(records[rec2])+1):
        l2 = records[rec2][j-1]
        s12 = match if l1 == l2 else lrmissmatch
        for k in range(1,len(records[rec3])+1):
            l3 = records[rec3][k-1]
            s13 = match if l1 == l3 else lrmissmatch
            s23 = match if l2 == l3 else srmissmatch
            gapxy = gap
            gapxz = gap
            gapyz = gap
            if k == len(records[rec3]):
                gapyz = 0
                gapxz = 0
            if j == len(records[rec2]):
                gapxy = 0
                gapyz = 0
            if i == len(records[rec1]):
                gapxy = 0
                gapxz = 0
            d = {'sx':mat[0][j][k], 'sy':mat[i][0][k], 'sz':mat[i][j][0], 'xyz': mat[i-1][j-1][k-1] + s12 + s23 + s13, 'xy': mat[i-1][j-1][k] + s12 + gapyz + gapxz, 'xz': mat[i-1][j][k-1] + s13 + gapxy + gapyz, 'yz': mat[i][j-1][k-1] + s23 + gapxy + gapxz, 'x': mat[i-1][j][k] + gapxy + gapxz, 'y': mat[i][j-1][k] + gapxy + gapyz, 'z': mat[i][j][k-1]+ gapxz + gapyz}
            #d = {'s':0, 'xyz':mat[i-1][j-1][k-1] + s12 + s23 + s13}
            a = max(d, key=d.get)
            traceback[i][j][k] = 0|bits[a];
            #print("ijk" + str(i) + ":"+ str(j) + ":" + str(k))
            #mat[i][j][k] = max(0, mat[i-1][j-1][k-1] + s12 + s23 + s13, mat[i-1][j-1][k] + s12 + gap,mat[i-1][j][k-1] + s13 + gap, mat[i][j-1][k-1] + s23 + gap, mat[i-1][j][k] + 2*gap, mat[i][j-1][k] + 2*gap, mat[i][j][k-1]+ 2*gap)
            mat[i][j][k] = d[a]
            #gap = gapo

#print(mat)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#surf(mat[1,:],mat[2,:],mat[3,1])

#print(traceback)
y = np.arange(0.0, len(records[rec1]))
x = np.arange(0.0, len(records[rec2]))
x,y = np.meshgrid(x,y)
#print(x)
#print(y)

fig = plt.figure()
ax = fig.gca(projection='3d')
#print(mat[0,1,:])
#print(mat[1,1,:])
#print(mat[1,:,:])
maxz = mat.max()
minz = mat.min()
#print(maxz)
#print(minz)
#for i in range(len(records[rec3])):
#    plt.subplot(len(records[rec3]), 1, i+1)
#    z = mat[i,:,:]
#    plt.imshow(z, cmap=plt.get_cmap('hot'), vmin= minz, vmax=maxz)
#plt.show()

print("mat max: "+ str(mat.max()))
#p.rint(np.argmax(mat,axis=1))
ind = np.unravel_index(np.argmax(mat, axis=None), mat.shape)
#print(ind)
#print(traceback[ind])
#print(mat[ind])
val = traceback[ind]
al = ["", "", ""]
#print(ind[2]-1)
#print(records[rec3][3])
i= 1
while val != 0 and val !=256 and val != 128:
    #print("step: " + str(i))
    i += 1
    #print(str(ind) + "\t" + str(bits.inv[traceback[ind]])+ "\t" + str(mat[ind[0]][ind[1]][ind[2]])+ str(records[rec1][ind[0]-1])) 
    if val & 64:
        al[0] = records[rec1][ind[0]-1] + al[0]
        al[1] = records[rec2][ind[1]-1] + al[1]
        al[2] = records[rec3][ind[2]-1] + al[2]
        ind = (ind[0]-1, ind[1]-1, ind[2]-1)
    elif val & 32:
        al[0] = "_" + al[0]
        al[1] = records[rec2][ind[1]-1] + al[1]
        al[2] = records[rec3][ind[2]-1] + al[2]
        ind = (ind[0], ind[1]-1, ind[2]-1)
    elif val & 16:
        al[0] = records[rec1][ind[0]-1] + al[0]
        al[1] = "_" + al[1]
        al[2] = records[rec3][ind[2]-1] + al[2]
        ind = (ind[0]-1, ind[1], ind[2]-1)
    elif val & 8:
        al[0] = records[rec1][ind[0]-1] + al[0]
        al[1] = records[rec2][ind[1]-1] + al[1]
        al[2] = "_" + al[2]
        ind = (ind[0]-1, ind[1]-1, ind[2])
    elif val & 4:
        al[0] = "_" + al[0]
        al[1] = "_" + al[1]
        al[2] = records[rec3][ind[2]-1] + al[2]
        ind = (ind[0], ind[1], ind[2]-1)
    elif val & 2:
        al[0] = "_" + al[0]
        al[1] = records[rec2][ind[1]-1] + al[1]
        al[2] = "_" + al[2]
        ind = (ind[0], ind[1]-1, ind[2])
    elif val & 1:
        al[0] = records[rec1][ind[0]-1] + al[0]
        al[1] = "_" + al[1]
        al[2] = "_" + al[2]
        ind = (ind[0]-1, ind[1], ind[2])
    val = traceback[ind]

print(al[1])
print(al[0])
print(al[2])
