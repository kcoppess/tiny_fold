import numpy as np

ms2_hairpin = 'ACAUGAGGAUCACCCAUGU'
ms2_hairpin_basepairs = np.array([[0,18], [1, 17], [2, 16], [3, 15], [4, 14], [6, 13], [7, 12]])

f = open('R101.txt','r')

labels = f.readline().split('\t')
NOC1 = labels.index('NumberOfClusters1')
S = labels.index('Sequence')
KDnl = labels.index('KDnoligand')

NumberOfClusters1 = []
Sequence = []
KDnoligand = []

b = f.readline().split('\t')
while b and len(b) > 1:
    if int(b[NOC1]) > 20 and (b[S].find(ms2_hairpin) != -1): # condition: number of clusters > 20 and the ms2 hairpin must be present
        NumberOfClusters1.append( int(b[NOC1]) )
        Sequence.append( b[S] )
        KDnoligand.append( float(b[KDnl]) )
    b = f.readline().split('\t')
print 'done filtering data according to NumberOfClusters1 > 20'

closing_bp_indices = []

for i in range(len(Sequence)):
    h = Sequence[i].index(ms2_hairpin)
    closing_bp_indices.append(ms2_hairpin_basepairs[0] + [h, h])

f.close()

