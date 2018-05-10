import numpy as np
import tinyfold as tf
import scipy.optimize as so
import matplotlib.pyplot as plt

rna = []
energies = [5.69, 6., 4.09, -7.09]
actual_bpp = []

w = 0.01

f = open("sequences_train.txt", 'r')
b = f.readline()
i = 0
data = 0.
while b:
    rna.append(tf.RNA(b, False, energies, True, True))
    data += len(b)**2
    actual_bpp.append(np.array(rna[i].get_bpp_full()))
    actual_bpp[i][np.isinf(actual_bpp[i])] = 0.
    b = f.readline()
    i += 1
f.close()

#guess = np.array([5., 4., 4., -6.])
guess = np.arange(-20,-16)

def cost(param, i, j):
    #param = np.array([p, 6., 4.09, -7.09])
    l2 = 0.
    for mm in range(i, j):
        rna[mm].update_energy(list(param))
        log_bpp = np.array(rna[mm].get_bpp_full())
        log_bpp[np.isinf(log_bpp)] = 0.
        term = (actual_bpp[mm] - log_bpp)
        l2 += np.sum(term ** 2)
    prior = guess - param
    l2 += w * np.dot(prior, prior)
    return l2

def cost_gradient(param, i, j):
    l2_grad = np.zeros(4)
    for mm in range(i, j):
        rna[mm].update_energy(list(param))
        grad = np.array(rna[mm].get_log_bpp_gradient_full())
        grad[np.isnan(grad)] = 0.
        log_bpp = np.array(rna[mm].get_log_bpp_full())
        log_bpp[np.isinf(log_bpp)] = 0.
        term = (actual_bpp[mm] - log_bpp)
        for ll in range(4):
            l2_grad[ll] += np.sum(-2 * term * grad[:,:,ll])
    prior = guess - param
    l2_grad += -2 * w * prior
    return l2_grad

param_iterations = []
prior_updates = []
def check(x):
    param_iterations.append(cost(x,0,50))
    new = guess - x
    print w*np.dot(new, new)
    prior_updates.append(w * np.dot(new, new))
    print str(cost(x,0,50)) + ' ' + str(cost_gradient(x, 0,50)) + ' ' + str(x)
    return

par = np.arange(-10,-6)
#par = np.arange(1,5)

def r(param):
    rna[0].update_energy(list(param))
    return rna[0].get_log_partition()#bpp(0,6)

def r_grad(param):
    rna[0].update_energy(list(param))
    return rna[0].get_log_gradient()#bpp_gradient(0,6)

#print so.check_grad(cost, cost_gradient, par, 0, 10)
print so.check_grad(r, r_grad, par)
#print r_grad(par) - so.approx_fprime(par, r, 1e-11)
'''
t = 0
while t < 1:
    for p in range(10):
        optimization = so.minimize(cost, par, args=(p*5, (p+1)*5), method='BFGS', jac=cost_gradient, tol=1e-15, callback=check)
        par = optimization.x
    t += 1
optimization = so.minimize(cost, par, args=(0, 50), method='BFGS', jac=cost_gradient, tol=1e-15, callback=check)
par = optimization.x
'''
optimization = so.minimize(cost, par, args=(0,50), method='Nelder-Mead', callback=check)
print optimization

rmsd = 0.
N = 0.
for w in range(50):
    log_bpp = np.array(rna[w].get_bpp_full())
    log_bpp[np.isinf(log_bpp)] = 0.
    N += np.count_nonzero(log_bpp)
    residual = np.absolute((actual_bpp[w] - log_bpp)/actual_bpp[w])
    residual[np.isnan(residual)] = 0.
    rmsd += np.sum(residual)
print rmsd/N

n = len(prior_updates)
k = range(n)

#e = np.linspace(-10, 15, 100)
#for w in [1.0, 0.1, 0.01]:
#    loss = []
#    for f in e:
#        loss.append(cost(f, 0, 50, w))
#    plt.plot(e, loss)

plt.plot(k, param_iterations, label = 'Loss')
plt.plot(k, prior_updates, label = 'Prior')
plt.legend()
plt.xlabel('Iteration')
plt.show() 
