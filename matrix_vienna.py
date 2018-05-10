import numpy as np
from subprocess import call

seq = open("experimental/R101_Sequence.txt", 'r')

for ii in range(1000): #12089):
    SEQUENCE = seq.readline().rstrip()
    NN = len(SEQUENCE)
    BPP_MATRIX = np.zeros((NN, NN))

    bpp = np.loadtxt('vienna_bpp_data/{}.txt'.format(SEQUENCE))
    for entry in bpp:
        BPP_MATRIX[int(entry[0])-1, int(entry[1])-1] = entry[2]
    print SEQUENCE
    np.savetxt('vienna_bpp_data/{}_bpp_matrix.txt'.format(SEQUENCE), BPP_MATRIX)
    call(["rm","vienna_bpp_data/{}.txt".format(SEQUENCE)])

seq.close()
