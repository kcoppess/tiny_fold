import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt
import filtering as fil

rna = []

alpha = 1.

w = 1e-22
energies = np.arange(-50,-46)

'''BPP TRAINING DATA'''
num_training_examples = 100
actual_bpp = [] # storing bpp for closing stem of hairpin for each input sequence
#energy_param = p.energies
sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
for i in range(num_training_examples):
    rna.append(tf.RNA(sequences[i], False, list(energies), True, True))
    actual_bpp.append((1e-9) / fil.KDnoligand[i])
print 'finished gathering training data'

guess = np.array([5.69, 6., 4.09, -7.09])

def cost(param, i, j):
    l2 = 0.
    for mm in range(i, j):
        rna[mm].update_energy(list(param))
        bp = fil.closing_bp_indices[mm]
        l2 += alpha * (actual_bpp[mm] - rna[mm].get_bpp(bp[0], bp[1])) ** 2
    prior = guess - param
    l2 += w * np.dot(prior, prior)
    return l2

def cost_gradient(param, i, j):
    l2_grad = np.zeros(4)
    for mm in range(i, j):
        rna[mm].update_energy(list(param))
        bp = fil.closing_bp_indices[mm]
        grad = np.array(rna[mm].get_bpp_gradient(bp[0], bp[1]))
        l2_grad += -2 * alpha * (actual_bpp[mm] - rna[mm].get_bpp(bp[0], bp[1])) * grad
    prior = guess - param
    l2_grad += -2 * w * prior
    return l2_grad

param_iterations = []
prior_updates = []
def check(x):
    param_iterations.append(cost(x,0,num_training_examples))
    new = guess - x
    prior_updates.append(w * np.dot(new, new))
    print str(cost(x,0,num_training_examples)) + ' ' + str(cost_gradient(x, 0, num_training_examples)) + ' ' + str(x)
    return


par = np.arange(-10,-6)
t = 0
while t < 1:
    for p in range(5):
        optimization = so.minimize(cost, par, args=(p*10, (p+1)*10), method='BFGS', jac=cost_gradient, tol=1e-20, callback=check)
        par = optimization.x
    t += 1
optimization = so.minimize(cost, par, args=(0, num_training_examples), method='BFGS', jac=cost_gradient, tol=1e-20, callback=check)
par = optimization.x

print optimization

final = []
for w in range(num_training_examples):
    bp = fil.closing_bp_indices[w]
    final.append(rna[w].get_bpp(bp[0],bp[1]))

y = np.array(final)
x = np.array(actual_bpp)

yhat = x
ybar = np.sum(y)/num_training_examples
ssreg = np.sum((yhat-ybar)**2)
sstot = np.sum((y - ybar)**2)
rsquared = 1 - (ssreg/sstot)
print rsquared

n = len(prior_updates)
k = range(n)

plt.plot(k, param_iterations, label = 'Loss')
plt.plot(k, prior_updates, label = 'Prior')
plt.legend()
plt.xlabel('Iteration')
plt.show()
