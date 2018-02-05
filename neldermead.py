import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

rna = []

alpha = 1.

w = 0.01
energies = np.arange(-10,-6)

'''BPP TRAINING DATA'''
num_training_examples = 100
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence
#energy_param = p.energies
sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
for i in range(num_training_examples):
    rna.append(tf.RNA(sequences[i], False, list(energies), True, False))
    actual_bpp.append(np.log(1e-9) - np.log(fil.KDnoligand[i]))
print 'finished gathering training data'

guess = np.array([5.69, 6., 4.09, -7.09])
#guess = np.array([0., 0., 0., 0.])

def cost(param, i, j):
    l2 = 0.
    for mm in range(i, j):
        rna[mm].update_energy(list(param))
        bp = fil.closing_bp_indices[mm]
        l2 += alpha * (actual_bpp[mm] - rna[mm].get_log_bpp(bp[0], bp[1])) ** 2
    prior = guess - param
    l2 += w * np.dot(prior, prior)
    return l2

param_iterations = []
prior_updates = []
def check(x):
    param_iterations.append(cost(x,0,num_training_examples))
    new = guess - x
    prior_updates.append(w * np.dot(new, new))
    print str(cost(x,0,num_training_examples)) + ' ' + str(x)
    return


par = np.arange(1,5)
optimization = so.minimize(cost, par, args=(0,num_training_examples), method='Nelder-Mead', callback=check)
print optimization

n = len(prior_updates)
k = range(n)

plt.plot(k, param_iterations, label = 'Loss')
plt.plot(k, prior_updates, label = 'Prior')
plt.legend()
plt.xlabel('Iteration')
plt.show()
