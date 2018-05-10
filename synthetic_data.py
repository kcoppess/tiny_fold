import tinyfold as tf
import numpy as np
import matplotlib.pyplot as plt

energies = [-0.93, -1.1, -1.33, -2.08, -2.11, -2.24, -2.35, -2.36, -3.26, -3.42, 4.09, 0.45]

Sequence = np.loadtxt("experimental/R101_Sequence.txt", dtype = 'string')
closing_bp_indices = np.loadtxt("experimental/R101_closing_bp_indices.txt", dtype = 'int')
ms2_hairpin_basepairs = np.loadtxt("experimental/R101_ms2_hairpin_basepairs.txt", dtype = 'int')

closingBP_kd = []
fullHP_kd = []
n = len(Sequence)

for mm in range(n):
    bp = closing_bp_indices[mm]
    closingBP_kd.append(10./tf.RNA(Sequence[mm], False, energies, True, False).get_bpp(bp[0],bp[1]))
    hp = ms2_hairpin_basepairs[mm]
    prob = 1.
    for ii in range(7):
        prob *= tf.RNA(Sequence[mm], False, energies, True, False).get_bpp(hp[2*ii], hp[2*ii + 1])
    fullHP_kd.append(10./prob)
    print mm
    #kd.append(np.log(1e-9) - tf.RNA(fil.Sequence[i], False, energies, True, False).get_log_bpp(bp[0],bp[1]))
    #log_exp_kd.append(np.log(fil.KDnoligand[i]))
np.savetxt("synthetic/R101_closing_basepair_KD.txt", closingBP_kd, delimiter='\t')
np.savetxt("synthetic/R101_full_ms2_KD.txt", fullHP_kd, delimiter='\t', fmt="%d")

#histo, edges = np.histogram(log_exp_kd + kd, bins='auto')

#plt.hist(log_exp_kd, bins=edges, label='Experimental')
#plt.hist(kd, bins=edges, label='Synthetic')
#plt.legend()
#plt.xlabel('log(Kd)')
#plt.ylabel('Frequency')
#plt.show()
