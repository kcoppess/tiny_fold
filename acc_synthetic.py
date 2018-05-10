import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt

rna = []
rna1 = []
rna2 = []

Sequence = np.loadtxt("experimental/R101_Sequence.txt", dtype = 'string', delimiter = '\t')
closing_bp_indices = np.loadtxt("experimental/R101_closing_bp_indices.txt", dtype = 'int', delimiter = '\t')
KDnoligand = np.loadtxt("synthetic/R101_closing_basepair_KD.txt", dtype = 'float', delimiter='\t')
#KDnoligand = np.loadtxt("synthetic/R101_full_ms2_KD.txt", dtype = 'float', delimiter='\t')
ms2_hairpin_basepairs = np.loadtxt("experimental/R101_ms2_hairpin_basepairs.txt", dtype = 'int', delimiter = '\t')

energies = np.array([-1.25618463, -0.39114066, -0.06533167, -3.57408081, -5.37062082,
       -4.10985842, -5.02943074, -3.55104567, -6.97317086, -8.67083103,
       12.12930106,  0.66978466])

'''BPP TRAINING DATA'''
num_training_examples = 100
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence

sequences = Sequence[:num_training_examples]#fil.Sequence[:num_training_examples]
for i in range(num_training_examples):
    #ms2_prob = 1.
    #bp = ms2_hairpin_basepairs[i]
    #for mm in range(7):
    #    ms2_prob *= tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(bp[2*mm], bp[2*mm + 1])
    bp = closing_bp_indices[i] #fil.closing_bp_indices[i]
    rna.append(np.log(10.) - tf.RNA(sequences[i], False, list(energies), True, False).get_log_bpp(bp[0], bp[1]))
    #rna.append(10./tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(bp[0], bp[1]))
    #rna.append(np.log(10.) - np.log(ms2_prob))
    actual_bpp.append(np.log(KDnoligand[i]))
    #actual_bpp.append(KDnoligand[i])
print 'finished gathering training data'

rna = np.array(rna)
actual_bpp = np.array(actual_bpp)

RMSD = np.mean((rna - actual_bpp)**2)

plt.title('R101 Synthetic w/ imperfect prior, hella far initial start')
plt.text(4, 17.5, "RMSD = %.2f"%RMSD)
#plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, label = '\"Uniform\" Kd', alpha=0.6, edgecolor = '')
plt.plot(actual_bpp, actual_bpp)
#plt.legend()
plt.xlabel('Synthetic Experimental log(Kd)')
plt.ylabel('Predicted log(Kd)')
plt.show()
