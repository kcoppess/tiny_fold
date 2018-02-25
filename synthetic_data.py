import filtering as fil
import tinyfold as tf
import numpy as np
import matplotlib.pyplot as plt

energies = [ 12.,   0.4,   0.3,  -1.]

kd = []
log_exp_kd = []
n = len(fil.Sequence)

for i in range(n):
    bp = fil.closing_bp_indices[i]
    kd.append(1e-9/tf.RNA(fil.Sequence[i], False, energies, True, False).get_bpp(bp[0],bp[1]))
    #kd.append(np.log(1e-9) - tf.RNA(fil.Sequence[i], False, energies, True, False).get_log_bpp(bp[0],bp[1]))
    #log_exp_kd.append(np.log(fil.KDnoligand[i]))
np.savetxt("synthetic_KD.txt", kd, delimiter='\t')
np.savetxt("eterna_sequences.txt", fil.Sequence, delimiter='\t', fmt="%s")
np.savetxt("ms2_hairpin_closing.txt", fil.closing_bp_indices, delimiter='\t', fmt="%d")

#histo, edges = np.histogram(log_exp_kd + kd, bins='auto')

#plt.hist(log_exp_kd, bins=edges, label='Experimental')
#plt.hist(kd, bins=edges, label='Synthetic')
#plt.legend()
#plt.xlabel('log(Kd)')
#plt.ylabel('Frequency')
#plt.show()
