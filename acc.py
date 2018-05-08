import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

Sequence = np.loadtxt("experimental/R101_Sequence.txt", dtype = 'string', delimiter = '\t')
KDnoligand = np.loadtxt("experimental/R101_KDnoligand.txt", dtype = 'float', delimiter = '\t')
closing_bp_indices = np.loadtxt("experimental/R101_closing_bp_indices.txt", dtype = 'int', delimiter = '\t')
ms2_hairpin_basepairs = np.loadtxt("experimental/R101_ms2_hairpin_basepairs.txt", dtype = 'int', delimiter = '\t')
vienna_ms2_bpp = np.loadtxt("R101_vienna_ms2_bpp.txt", dtype = 'float')

rna = []
rna1 = []
rna2 = []

#energies = np.array([ 13.53115792,   0.41636607,   4.45098914,  -0.63286917]) # trained on uniform filtered data, w = 0.01
#energies2 = np.array([ 13.58526261,   3.47537565,   2.04322166,  -2.29645279]) # w = 0.01

energies = np.array([ -2.62797039,  -0.94625975,   1.55446311,   0.48276247,
        -5.44527292,  -6.11612285,  -0.63093985,  -0.11027562,
       -12.61154124,   2.9198702 ,   1.2749916 ,  -3.2458907 ])
#energies = np.array([-0.93, -1.1, -1.33, -2.08, -2.11, -2.24, -2.35, -2.36, -3.26, -3.42, 4.09, 0.45])

'''BPP TRAINING DATA'''
num_training_examples = 100
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence

trouble_sequences = []

sequences = Sequence[:num_training_examples]#fil.Sequence[:num_training_examples]
for mm in range(num_training_examples):
    #bp = closing_bp_indices[i] #fil.closing_bp_indices[i]
    #rna.append(np.log(1e-9) - tf.RNA(sequences[i], False, list(energies), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    ms2_prob = 1.
    bp = ms2_hairpin_basepairs[mm]
    for ii in range(7):
        ms2_prob *= tf.RNA(sequences[mm], False, list(energies), True, False).get_bpp(bp[2*ii], bp[2*ii + 1])
    rna.append(ms2_prob)
    if ((ms2_prob < 1e-5) and (vienna_ms2_bpp[mm] > 1e-5)):
        trouble_sequences.append(sequences[mm])
    actual_bpp.append(vienna_ms2_bpp[mm])
    #rna.append(np.log(ms2_prob))
    #actual_bpp.append(np.log(vienna_ms2_bpp[mm]))
    #rna.append(ms2_prob)
    #actual_bpp.append(10./KDnoligand[mm])
print 'finished gathering training data'

#np.savetxt('trouble.txt', trouble_sequences, fmt = '%s')
#np.savetxt('guess_kd.txt', rna)

plt.title('Full MS2 Hairpin Probability')
#plt.title('R101 Training Set (BFGS + 2 runs of NM)')
#plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, label = '\"Uniform\" Kd', alpha = 0.6, edgecolor = '')
plt.plot(actual_bpp, actual_bpp)
#plt.legend()
#plt.xlabel('Vienna')
#plt.ylabel('tinyFold')
plt.xlabel('Experimental log(Kd)')
plt.ylabel('Predicted log(Kd)')
plt.show()

