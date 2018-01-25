import numpy as np
import tinyfold as tf
import scipy.optimize as so

rna = []
energies = [5.69, 6., 4.09, -7.09]
actual_bpp = []

f = open("sequences_train.txt", 'r')
b = f.readline()
i = 0
while b:
    rna.append(tf.RNA(b, False, energies, True, True))
    actual_bpp.append(np.array(rna[i].get_bpp_full()))
    b = f.readline()
    i += 1
f.close()

def cost(param):
    l2 = 0.
    for mm in range(50):
        rna[mm].update_energy(param)
        l2 += np.sum((actual_bpp[mm] - np.array(rna[mm].get_bpp_full())) ** 2)
    return l2

def cost_gradient(param):
    l2_grad = np.zeros(4)
    for mm in range(50):
        rna[mm].update_energy(param)
        l2_grad += np.sum(np.sum(- 2 * np.dot((actual_bpp[mm] - np.array(rna[mm].get_bpp_full())), np.array(rna[mm].get_bpp_gradient_full())), axis=0), axis = 0)
    return l2_grad

print cost_gradient(energies)
