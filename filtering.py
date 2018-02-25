import numpy as np
import matplotlib.pyplot as plt

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
    if int(b[NOC1]) > 20 and (b[S].find(ms2_hairpin) != -1): #and (float(b[KDnl]) < 25.): # condition: number of clusters > 20 and the ms2 hairpin must be present
        NumberOfClusters1.append( int(b[NOC1]) )
        Sequence.append( b[S] )
        KDnoligand.append( float(b[KDnl]) )
    b = f.readline().split('\t')
print 'done filtering data according to NumberOfClusters1 > 20'

f.close()
'''
hist, edges = np.histogram(KDnoligand, bins='auto')

print np.count_nonzero(hist)
print len(hist)

NumberOfClusters11= []
Sequence1 = []
KDnoligand1 = []

binnum = len(edges)
bn = 1

while bn < binnum:
    for i in range(len(KDnoligand)):
        if (KDnoligand[i] <= edges[bn]) and (KDnoligand[i] > edges[bn-1]):
            KDnoligand1.append(KDnoligand[i])
            Sequence1.append(Sequence[i])
            NumberOfClusters11.append(NumberOfClusters1[i])
            break
    bn += 1


KDnoligand = KDnoligand1
Sequence = Sequence1
NumberOfClusters1 = NumberOfClusters11
print len(KDnoligand)
'''

closing_bp_indices = []
for i in range(len(Sequence)):
    h = Sequence[i].index(ms2_hairpin)
    closing_bp_indices.append(ms2_hairpin_basepairs[0] + [h, h])
'''
plt.hist(KDnoligand, bins='auto')
plt.title('Experimental Kd Distribution')
plt.xlabel('Kd')
plt.ylabel('Frequency')
plt.show()
'''
