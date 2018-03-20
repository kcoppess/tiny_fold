import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

rna = []
rna1 = []
rna2 = []

#energies = np.array([ 13.53115792,   0.41636607,   4.45098914,  -0.63286917]) # trained on uniform filtered data, w = 0.01
#energies2 = np.array([ 13.58526261,   3.47537565,   2.04322166,  -2.29645279]) # w = 0.01

energies = np.array([ 1.24716782,  6.68217972, -2.30364506, -1.68943857,  6.6996317 ,
        6.35152126,  4.96016927,  6.22891496,  4.33263617,  2.66201513,
        8.59493876,  4.42860694])

'''BPP TRAINING DATA'''
num_training_examples = 30
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence

sequences = fil.Sequence[:num_training_examples]#fil.Sequence[:num_training_examples]
for i in range(num_training_examples):
    bp = fil.closing_bp_indices[i] #fil.closing_bp_indices[i]
    rna.append(np.log(1e-9) - tf.RNA(sequences[i], False, list(energies), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    #rna1.append(tf.RNA(sequences[i], False, list(energies1), True, False).get_log_bpp(bp[0], bp[1]))
    #rna2.append(np.log(1e-9) - tf.RNA(sequences[i], False, list(energies2), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(int(bp[0]), int(bp[1])))
    #rna1.append(1e-9/tf.RNA(sequences[i], False, list(energies1), True, False).get_bpp(bp[0], bp[1]))
    #rna2.append(1e-9/tf.RNA(sequences[i], False, list(energies2), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(fil.KDnoligand[i])
    actual_bpp.append(np.log(fil.KDnoligand[i]))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(fil.KDnoligand[i])
print 'finished gathering training data'

plt.title('R101 Training Set')
#plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, label = '\"Uniform\" Kd')
plt.plot(actual_bpp, actual_bpp)
#plt.legend()
plt.xlabel('Experimental log(Kd)')
plt.ylabel('Predicted log(Kd)')
plt.show()
