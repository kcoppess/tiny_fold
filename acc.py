import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

rna = []

energies = np.array([ 13.16556948,   2.48849498,   1.60773765,  -1.74590936])
#energies = np.array([ 12.89831853,  11.62084721,   1.36910553,  -1.25232382]) #50
#energies = np.array([ 13.31310791,   6.13255856,   5.69664088,  -5.33310487]) #50 larger w

'''BPP TRAINING DATA'''
num_training_examples = 1000
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence
#energy_param = p.energies
sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
for i in range(num_training_examples):
    bp = fil.closing_bp_indices[i]
    rna.append(tf.RNA(sequences[i], False, list(energies), True, False).get_log_bpp(bp[0], bp[1]))
    actual_bpp.append(np.log(1e-9) - np.log(fil.KDnoligand[i]))
    #rna.append(1e-9/tf.RNA(sequences[i], False, list(energies), True, False).get_bpp(bp[0], bp[1]))
    #actual_bpp.append(fil.KDnoligand[i])
print 'finished gathering training data'

plt.title('Training Set')
plt.scatter(actual_bpp, rna)
plt.plot(rna, rna)
plt.legend()
plt.xlabel('Experiment')
plt.show()
