import numpy as np
import matplotlib.pyplot as plt
import parameters as p
import base_pair_probabilities as bpp
import g_matrix as gm
import filtering as fil
 
num_training_examples = 50
training_data = [] # storing bpp for closing stem of hairpin for each input sequence
g_loop = p.g_loop
sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
for i in range(num_training_examples):
    training_data.append((1e-9) / fil.KDnoligand[i])
print 'finished gathering training data'

# plotting for just 2 params
def cost(param): # param = [g_AU, g_stack]
    sum_sq = 0.

    for i in range(len(training_data)):
        M = len(sequences[i])
        g_base_pair = gm.generator(sequences[i], param[0], p.energies[1], p.energies[2], M)
    
        bp = fil.closing_bp_indices[i]
        bpp_matrix = bpp.mccaskill_linear(g_base_pair, g_loop, param[1], M)
        residual = training_data[i] - bpp_matrix[bp[0], bp[1]]
        sum_sq += residual**2
    return sum_sq

N = 100
g_stack = np.linspace(-15, 20, N)
g_AU = np.linspace(-15, 20, N)
cost_values = np.zeros((N,N))

for i in range(N):
    for j in range(N):
        print i, j
        params = np.array([g_AU[i], g_stack[j]])
        cost_values[i,j] = cost(params)


print 'g_stack = '+str(g_stack)
print '------------------------------'
print '------------------------------'
print '------------------------------'
print '------------------------------'
print 'g_AU = '+str(g_AU)
print '------------------------------'
print '------------------------------'
print '------------------------------'
print '------------------------------'
print cost_values

plt.imshow(cost_values, extent=[min(g_stack), max(g_stack), min(g_AU), max(g_AU)], origin='lower')
plt.xlabel('g_stack')
plt.ylabel('g_AU')
plt.colorbar()
plt.show()

