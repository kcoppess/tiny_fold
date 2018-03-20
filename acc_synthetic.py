import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

rna = []
rna1 = []
rna2 = []

Sequence = np.loadtxt("eterna_sequences.txt", dtype = 'string', delimiter='\t')
KDnoligand = np.loadtxt("synthetic_kd.txt", dtype = 'float', delimiter='\t')
closing_bp_indices = np.loadtxt("ms2_hairpin_closing.txt", delimiter='\t')

energies = np.array([ 1.59302437,  3.60824333,  5.42886472,  7.57214443,  8.1729206 ,
        6.09160341,  7.73909905,  8.37004938, 13.75500795,  6.87480214,
        1.23732032,  0.95195283])
#energies = np.array([ 11.99979688,   0.39980603,   0.29989286,  -0.99984654]) #synthetic data
#energies1 = np.array([ 12.58417779,   2.94746754,   6.77340698,   0.83004456]) # trained on even more filtered data, w = 0.01
#energies2 = np.array([ 13.58526261,   3.47537565,   2.04322166,  -2.29645279]) # w = 0.01

'''BPP TRAINING DATA'''
num_training_examples = 30
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence

sequences = Sequence[:num_training_examples]#fil.Sequence[:num_training_examples]
for i in range(num_training_examples):
    bp = closing_bp_indices[i] #fil.closing_bp_indices[i]
    rna.append(np.log(1e-9) - tf.RNA(sequences[i], False, list(energies), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    #rna1.append(tf.RNA(sequences[i], False, list(energies1), True, False).get_log_bpp(bp[0], bp[1]))
    #rna2.append(tf.RNA(sequences[i], False, list(energies2), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(int(bp[0]), int(bp[1])))
    #rna1.append(1e-9/tf.RNA(sequences[i], False, list(energies1), True, False).get_bpp(bp[0], bp[1]))
    #rna2.append(1e-9/tf.RNA(sequences[i], False, list(energies2), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(KDnoligand[i])
    actual_bpp.append(np.log(KDnoligand[i]))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(fil.KDnoligand[i])
print 'finished gathering training data'

plt.title('R101 Synthetic w/ perfect prior')
#plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, label = '\"Uniform\" Kd')
plt.plot(actual_bpp, actual_bpp)
#plt.legend()
plt.xlabel('Synthetic Experimental log(Kd)')
plt.ylabel('Predicted log(Kd)')
plt.show()
