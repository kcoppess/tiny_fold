import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil
import numpy.random as rd

rna = []

w = 1e-2
energies = np.arange(-10,-6)

#Sequence = fil.Sequence
#KDnoligand = fil.KDnoligand
#closing_bp_indices = fil.closing_bp_indices
Sequence = np.loadtxt("eterna_sequences.txt", dtype = 'string', delimiter='\t')
KDnoligand = np.loadtxt("synthetic_kd.txt", dtype = 'float', delimiter='\t')
closing_bp_indices = np.loadtxt("ms2_hairpin_closing.txt", delimiter='\t')

'''BPP TRAINING DATA'''
num_training_examples = 30
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence
#energy_param = p.energies
sequences = Sequence[:num_training_examples] #fil.Sequence[:num_training_examples]
for i in range(num_training_examples):
    rna.append(tf.RNA(sequences[i], False, list(energies), True, False))
    #actual_bpp.append(1e-9/fil.KDnoligand[i])
    #actual_bpp.append(np.log(1e-9) - np.log(fil.KDnoligand[i]))
    #actual_bpp.append(1e-9/KDnoligand[i])
    actual_bpp.append(np.log(1e-9) - np.log(KDnoligand[i]))
print 'finished gathering training data'

#guess = np.array([5.69, 6., 4.09, -7.09])
#guess = np.array([5., 5., 3.5, -6.])
guess = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])

def cost(param, indices, i):
    l2 = 0.
    for mm in indices: #range(i, j):
        rna[mm].update_energy(list(param))
        bp = closing_bp_indices[mm] #fil.closing_bp_indices[mm]
        l2 += (actual_bpp[mm] - rna[mm].get_log_bpp(int(bp[0]), int(bp[1]))) ** 2
    prior = guess - param
    l2 += w * np.dot(prior, prior)
    return l2

param_iterations = []
prior_updates = []
def check(x):
    param_iterations.append(cost(x,np.arange(num_training_examples),0))
    new = guess - x
    prior_updates.append(w * np.dot(new, new))
    #print w*np.dot(new,new)
    print str(cost(x,np.arange(num_training_examples),0)) + ' ' + str(x)
    return

#error = []
#for t in range(5, num_training_examples, 5):
#print 't = '+str(t)+ ' --------------------------------------'
err = 0.
iterations = 1
for k in range(iterations):
    #print k
    ind = np.arange(num_training_examples)
    #ind = rd.randint(num_training_examples, size = t)
    par = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.]
    #par = range(1, 13)
    optimization = so.minimize(cost, par, args=(ind, 0), method='Nelder-Mead', callback=check, options={'maxfev': 10000})
    print optimization
    
    final = []
    subset_actual = []
    for w in ind:
        bp = closing_bp_indices[w] #fil.closing_bp_indices[w]
        final.append(rna[w].get_log_bpp(int(bp[0]), int(bp[1])))
        subset_actual.append(actual_bpp[w])
    
    y = np.array(final)
    x = np.array(subset_actual)
    
    err += np.mean(np.absolute((y - x)/x))/iterations
    #error.append(err)
print err

#np.savetxt("error.txt", error, delimiter=',')
#yhat = x
#ybar = np.sum(y)/num_training_examples
#ssreg = np.sum((yhat-ybar)**2)
#sstot = np.sum((y - ybar)**2)
#rsquared = (ssreg/sstot)
#print rsquared

n = len(prior_updates)
k = range(n)

plt.plot(k, param_iterations, label = 'Loss')
plt.plot(k, prior_updates, label = 'Prior')
plt.legend()
plt.xlabel('Iteration')
plt.show()
