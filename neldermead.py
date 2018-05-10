import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
#import filtering as fil
import numpy.random as rd

rna = []

w = 1e-3
energies = np.arange(-10,2)

closing_bp_indices = np.loadtxt("experimental/R101_closing_bp_indices.txt", dtype = 'int', delimiter = '\t')
ms2_hairpin_basepairs = np.loadtxt("experimental/R101_ms2_hairpin_basepairs.txt", dtype = 'int', delimiter = '\t')
vienna_ms2_bpp = np.loadtxt("R101_vienna_ms2_bpp.txt", dtype = 'float')
vienna_closing_bpp = np.loadtxt("R101_vienna_closing_bpp.txt", dtype = 'float')
Sequence = np.loadtxt("experimental/R101_Sequence.txt", dtype = 'string', delimiter = '\t')
KDnoligand = np.loadtxt("experimental/R101_KDnoligand.txt", dtype = 'float', delimiter = '\t')
#KDnoligand = np.loadtxt("synthetic/R101_closing_basepair_KD.txt", dtype = 'float', delimiter='\t')

'''BPP TRAINING DATA'''
num_training_examples = 100
actual_bpp = [] # storing log(bpp) for closing stem of hairpin for each input sequence
#energy_param = p.energies
sequences = Sequence[:num_training_examples] #fil.Sequence[:num_training_examples]
for i in range(num_training_examples):
    rna.append(tf.RNA(sequences[i], False, list(energies), True, False))
    actual_bpp.append(10./KDnoligand[i])
    #actual_bpp.append(np.log(10.) - np.log(KDnoligand[i])) # FIXME
    #actual_bpp.append(np.log(vienna_ms2_bpp[i])) # FIXME
    #actual_bpp.append(vienna_closing_bpp[i])
    #actual_bpp.append(vienna_ms2_bpp[i]) # FIXME
print 'finished gathering training data'


#guess = np.array([5.69, 6., 4.09, -7.09])
#guess = np.array([5., 5., 3.5, -6.])
#guess = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
#guess = np.zeros(12)
#random_vec = np.random.uniform(-1, 1, size=12)
#offset = 0.5*random_vec / np.linalg.norm(random_vec)
#guess = offset + np.array([-0.93, -1.1, -1.33, -2.08, -2.11, -2.24, -2.35, -2.36, -3.26, -3.42, 4.09, 0.45])
guess = np.array([-0.93, -1.1, -1.33, -2.08, -2.11, -2.24, -2.35, -2.36, -3.26, -3.42, 4.09, 0.45])

def cost(param, indices, i):
    l2 = 0.
    for mm in indices: #range(i, j):
        rna[mm].update_energy(list(param))
        '''
        ms2_prob = 1.
        bp = ms2_hairpin_basepairs[mm]
        for ii in range(7):
            #print rna[mm].get_bpp(bp[2*ii], bp[2*ii + 1])
            ms2_prob *= rna[mm].get_bpp(bp[2*ii], bp[2*ii + 1])
            #print ":: " + str(rna[mm].get_bpp(bp[2*ii], bp[2*ii + 1]))
        #print np.log(ms2_prob)
        #print ":: " + str(actual_bpp[mm])
        l2 += (actual_bpp[mm] - ms2_prob) ** 2
        #l2 += (actual_bpp[mm] - np.log(ms2_prob)) ** 2
        #print ms2_prob
        '''
        bp = closing_bp_indices[mm] #fil.closing_bp_indices[mm]
        l2 += (actual_bpp[mm] - rna[mm].get_bpp(bp[0], bp[1])) ** 2
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
    initial = np.array([-0.93, -1.1, -1.33, -2.08, -2.11, -2.24, -2.35, -2.36, -3.26, -3.42, 4.09, 0.45])
    random_vec = np.random.uniform(-1, 1, size=len(initial))
    offset = 0.3*random_vec / np.linalg.norm(random_vec)
    par = initial + offset
    #par = np.array([-1.94370382, -0.13398913, -2.22103728, -0.92180998, -2.2727062 ,
    #   -1.19918207,  0.20761484, -1.59890573, -5.95664947, -0.3693114 ,
    #   -5.49697268, -1.28644052])
    optimization = so.minimize(cost, par, args=(ind, 0), method='Nelder-Mead', callback=check, options={'maxfev': 100000})
    print optimization

n = len(prior_updates)
k = range(n)

print initial + offset

print "experimental data, closing bpp"
'''
plt.plot(k, param_iterations, label = 'Loss')
plt.plot(k, prior_updates, label = 'Prior')
#plt.legend()
plt.xlabel('Iteration')
plt.show()
'''
