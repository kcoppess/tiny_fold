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
vienna_closing_bpp = np.loadtxt("R101_vienna_closing_bpp.txt", dtype = 'float')

rna = []
rna1 = []
rna2 = []

'''VIENNA DATA'''
# w = 1e-3
energies = np.array([-0.2427841 , -4.03141534, -1.09296522, -0.77945274, -1.45638276,
       -4.14362266,  0.28071853, -7.08002598, -1.19004444, -0.59074734,
        9.00680994,  0.00927074])
energies = np.array([-1.03818521e+00, -2.61586621e+00, -1.24002835e+00, -4.35664370e-02,
       -1.48861769e+00, -2.91078949e+00,  1.85821572e-01, -5.31165743e+00,
       -1.76672998e+00,  6.07027829e-01,  6.27049861e+00,  3.95830197e-03])
energies = np.array([-2.48962188e-01, -5.39598613e+00, -1.12410538e+00, -4.96916092e+00,
       -2.11125238e+00, -5.24480162e+00,  3.05836003e+00,  4.53132066e+00,
       -3.40350145e+00, -9.86159175e-02,  1.26395945e+01,  8.43958355e-04])
energies = np.array([-1.44822386e+00, -1.37366663e+00,  8.46557630e-02,  4.63638017e-01,
       -4.07925758e+00, -3.85734843e+00, -1.38211145e+00, -1.14268096e+01,
       -2.22046341e+00, -7.36998587e-01,  1.13195017e+01, -8.83318530e-03])
energies = np.array([-1.96782894e+00, -5.81045954e+00, -1.94943873e+00, -2.57391133e+00,
       -2.36506469e+00, -5.90909851e+00,  5.56624666e-02,  2.16528320e+00,
       -3.83135553e+00, -2.30475908e-01,  1.41680852e+01,  1.10841815e-02])
'''R101 DATA'''
energies = np.array([])
#energies = np.array([-0.93, -1.1, -1.33, -2.08, -2.11, -2.24, -2.35, -2.36, -3.26, -3.42, 4.09, 0.45])

'''BPP TRAINING DATA'''
num_training_examples = 100
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence

trouble_sequences = []

sequences = Sequence[:num_training_examples]#fil.Sequence[:num_training_examples]
for mm in range(num_training_examples):
    bp = closing_bp_indices[mm] #fil.closing_bp_indices[i]
    rna.append(tf.RNA(sequences[mm], False, list(energies), True, False).get_bpp(bp[0], bp[1]))
    actual_bpp.append(vienna_closing_bpp[mm])
    #rna.append(np.log(10.) - tf.RNA(sequences[mm], False, list(energies), True, False).get_log_bpp(bp[0], bp[1]))
    #actual_bpp.append(np.log(vienna_closing_bpp[mm]))
    #ms2_prob = 1.
    #bp = ms2_hairpin_basepairs[mm]
    #for ii in range(7):
    #    ms2_prob *= tf.RNA(sequences[mm], False, list(energies), True, False).get_bpp(bp[2*ii], bp[2*ii + 1])
    #rna.append(ms2_prob)
    #if ((ms2_prob < 1e-5) and (vienna_ms2_bpp[mm] > 1e-5)):
    #    trouble_sequences.append(sequences[mm])
    #actual_bpp.append(vienna_ms2_bpp[mm])
    #rna.append(np.log(ms2_prob))
    #actual_bpp.append(np.log(vienna_ms2_bpp[mm]))
    #rna.append(ms2_prob)
    #actual_bpp.append(10./KDnoligand[mm])
print 'finished gathering training data'

rna = np.array(rna)
actual_bpp = np.array(actual_bpp)

RMSD = np.mean((rna - actual_bpp)**2)

#np.savetxt('trouble.txt', trouble_sequences, fmt = '%s')
#np.savetxt('guess_kd.txt', rna)

plt.text(0.1, 0.8, "RMSD = %.2f"%RMSD)
#plt.title('R101 Training Set (BFGS + 2 runs of NM)')
#plt.scatter(actual_bpp, rna2, c='m', label = 'Full')
#plt.scatter(actual_bpp, rna1, c='c', label = 'Kd < 25')
plt.scatter(actual_bpp, rna, alpha = 0.6, edgecolor = '')
plt.plot(actual_bpp, actual_bpp)
#plt.legend()
plt.xlabel('Vienna BPP')
plt.ylabel('tinyFold BPP')
#plt.xlabel('Vienna log(Kd)')
#plt.ylabel('Predicted log(Kd)')
plt.show()

