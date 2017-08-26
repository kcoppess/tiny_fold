# modified: August 26, 2017
# using Adam method
import numpy as np
import log_partition as z
import training_sequence_generator as tr 
import random

# returns gradient of residual squared for particular point
# NOTE new cost function: sum( (Qi - Q(xi,B))^2 ) + alpha*B^2
def gradient(energy_param,alpha,sequence,partition):
    N = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring '-'

    grad = 2*alpha*energy_param
    
    # ignoring steric effects, generating base-pair free energy matrix
    g_base_pair = np.zeros((N,N))
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = z.free_energy_pair(sequence[m], sequence[n], energy_param[0], energy_param[1], energy_param[2]) # bp1, bp2, g_AU, g_GU, g_GC

    if is_circular:
        pred_partition = z.circular_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N)
        for j in range(5):
            gr = z.circular_sequence_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, g
            grad[j] += (-2)*(partition - pred_partition) * gr 
    else:
        pred_partition = z.linear_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N)
        for j in range(5):
            gr = z.linear_sequence_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
            grad[j] += (-2)*(partition - pred_partition) * gr
    #grad[1] = grad[2] = grad[3] = grad[4] = 0.
    return grad

def convergence(param, prev_param): 
    if np.linalg.norm(param - prev_param) < 1e-10:
        return True
    else:
        return False

param = np.array([1.,2.,3.,4.,5.])
prev_param = np.zeros(5)
index = range(len(tr.training_data))

a = 0.0001

alpha = 0.0001
beta1 = 0.9
beta2 = 0.999
m = np.zeros(5)
v = np.zeros(5)
t = np.zeros(5)

k = 1

K = []
iteration_param = [] # list of updated parameters after each iteration

while not convergence(param, prev_param) and k < 4001:
    random.shuffle(index) # randomly shuffling indexes for each pass through the data
    print k
    for i in index:
        prev_param = np.array(param)
        t += 1
        grad = gradient(param, a, tr.sequences[i], tr.training_data[i])
        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad**2
        mhat = m / (1 - beta1**t)
        vhat = v / (1 - beta2**t)
        param += (-alpha) * mhat / (np.sqrt(vhat) + 1e-8)
        #print tr.sequences[i], tr.training_data[i]
        #print 'Gradient: '+str(grad)
        #print 'Prev:    '+str(prev_param)
        #print 'Updated: '+str(param)
        #print np.linalg.norm(param-prev_param)
        #print '-------------------'
    K.append(k)
    iteration_param.append(list(param))
    k += 1
print param
