import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

rna = []
rna1 = []
rna2 = []

#energies = np.array([ 5.68986683,  5.99983762,  4.08992768, -7.08986771]) #synthetic data
energies = np.array([ 13.53115792,   0.41636607,   4.45098914,  -0.63286917]) # trained on uniform filtered data, w = 0.01
#energies1 = np.array([ 12.58417779,   2.94746754,   6.77340698,   0.83004456]) # trained on even more filtered data, w = 0.01
##energies = np.array([ 12.69638911,   2.4916732 ,   6.7655894 ,   1.06377874]) # trained on filtered data, w = 0.01
##energies = np.array([ 12.69553677,   2.48808272,   7.92487499,   2.32077606]) # trained on filtered data, w = 0.001
##energies = np.array([ 13.58514293,   3.4746541 ,   2.04299684,  -2.2962257 ]) # w = 0.001
##energies = np.array([ 13.75490208,   5.49259234,   4.2242013 ,  -4.16502061]) # w = 0.1
energies2 = np.array([ 13.58526261,   3.47537565,   2.04322166,  -2.29645279]) # w = 0.01
##energies1 = np.array([ 13.25114318,   3.18594587,   1.70107373,  -1.77127592]) # 100 examples, w = 0.01

'''BPP TRAINING DATA'''
num_training_examples = 8000
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence

sequences = fil.Sequence[:num_training_examples]#fil.Sequence[:num_training_examples]
for i in range(num_training_examples):
    bp = fil.closing_bp_indices[i] #fil.closing_bp_indices[i]
    rna.append(np.log(1e-9) - tf.RNA(sequences[i], False, list(energies), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    #rna1.append(tf.RNA(sequences[i], False, list(energies1), True, False).get_log_bpp(bp[0], bp[1]))
    rna2.append(np.log(1e-9) - tf.RNA(sequences[i], False, list(energies2), True, False).get_log_bpp(int(bp[0]), int(bp[1])))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(int(bp[0]), int(bp[1])))
    #rna1.append(1e-9/tf.RNA(sequences[i], False, list(energies1), True, False).get_bpp(bp[0], bp[1]))
    #rna2.append(1e-9/tf.RNA(sequences[i], False, list(energies2), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(fil.KDnoligand[i])
    actual_bpp.append(np.log(fil.KDnoligand[i]))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(fil.KDnoligand[i])
print 'finished gathering training data'

plt.title('R101 Training Set')
plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, label = '\"Uniform\" Kd')
plt.plot(actual_bpp, actual_bpp)
plt.legend()
plt.xlabel('Experimental log(Kd)')
plt.ylabel('Predicted log(Kd)')
plt.show()
