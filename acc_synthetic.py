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

energies = np.array([-0.00961283,  0.10465754,  0.25906695,  0.22889634,  0.49542794,  0.59705241,
  0.68148017,  0.57226634,  0.82583153,  1.03634249,  1.10152011,  1.20060662]) # imperfect prior 8000 close initial
#energies = np.array([-0.00368309,  0.12300926,  0.27529714,  0.30213592,  0.4970385,   0.5505612,
#  0.69836815,  0.55733555,  0.92233248,  1.14318185,  1.10053463,  1.21615434]) # perfect prior 8000 close initial
#energies = np.array([-0.59826039, -1.97512995, -0.35859398, -0.56573092,  1.79617267,
#       -0.37297762,  0.81471193, -0.2519902 , -0.17115048, -0.73916539,
#        1.20170385,  1.61639218]) # perfect prior 8000
#energies = np.array([-0.18518726, -1.48659774, -1.41313968,  0.04836343,  0.52431478,
#       -2.56340403,  0.38445788, -0.5796633 , -0.1225757 ,  0.7479511 ,
#        1.17450923,  3.58602338]) # perfect prior
#energies = np.array([-0.62611138, -1.90739849, -0.30264324, -0.48317788,  1.70992573,
#       -0.3327077 ,  0.68667703, -0.31308698, -0.12954505, -0.7700525 ,
#        1.20060026,  1.5928439 ]) # imperfect prior 8000
#energies = np.array([-2.02587778, -1.6336197 , -0.47941516, -1.04056735,  3.50311667,
#       -0.97938649,  0.39213955, -1.08824477, -0.02494981,  1.64171906,
#        1.10184002,  2.11486854]) # imperfect prior
#energies = np.array([ 11.99979688,   0.39980603,   0.29989286,  -0.99984654]) #synthetic data
#energies1 = np.array([ 12.58417779,   2.94746754,   6.77340698,   0.83004456]) # trained on even more filtered data, w = 0.01
#energies2 = np.array([ 13.58526261,   3.47537565,   2.04322166,  -2.29645279]) # w = 0.01

'''BPP TRAINING DATA'''
num_training_examples = 8000
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

plt.title('R101 Synthetic w/ imperfect prior, close initial start')
#plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, label = '\"Uniform\" Kd')
plt.plot(actual_bpp, actual_bpp)
#plt.legend()
plt.xlabel('Synthetic Experimental log(Kd)')
plt.ylabel('Predicted log(Kd)')
plt.show()
