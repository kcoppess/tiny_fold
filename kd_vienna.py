import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

num = 130#12089

seq = open("experimental/R101_Sequence.txt", 'r')
cbp = np.loadtxt("experimental/R101_closing_bp_indices.txt")
mbp = np.loadtxt("experimental/R101_ms2_hairpin_basepairs.txt")
kd_experimental = np.loadtxt("experimental/R101_KDnoligand.txt")

kd = []

for ii in range(num):#12089):
    SEQUENCE = seq.readline().rstrip()
    CBP = cbp[ii,:]
    MBP = mbp[ii,:]
    
    matrix = np.loadtxt('vienna_bpp_data/{}_bpp_matrix.txt'.format(SEQUENCE))

    hairpin_prob = 1.#matrix[int(CBP[0]), int(CBP[1])]
    for mm in range(7):
        hairpin_prob *= matrix[int(MBP[2*mm]), int(MBP[2*mm + 1])]
    kd.append(hairpin_prob)
    #print hairpin_prob
    print ii

#np.savetxt("R101_vienna_ms2_bpp.txt", np.array(kd))

seq.close()

x = np.array(kd)
y = 1./np.array(kd_experimental[:num])

xy = np.vstack([np.array(kd), kd_experimental[:num]])
z = gaussian_kde(xy)(xy)

idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

print "plotting..."

#plt.scatter(x, y, c=z, cmap='jet')
plt.scatter(x, y, alpha=0.6, edgecolor = '')
plt.xlabel("Full Hairpin Probability (Vienna)")
plt.ylabel("Experimental 1/KD")
plt.show()

