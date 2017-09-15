# modified: Septemeber 14, 2017
# using Adam method
import numpy as np
import log_partition as log_z
#import random_sequence_generator as rsg 
import parameters as p
import random

# returns gradient of residual squared for particular point
# NOTE new cost function: sum( (Qi - Q(xi,B))^2 ) + alpha*B^2
def gradient(energy_param, sequence, partition):
    N = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring '-'

    grad = 2 * p.w * energy_param
    
    # ignoring steric effects, generating base-pair free energy matrix
    g_base_pair = np.zeros((N,N))
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = log_z.free_energy_pair(sequence[m], sequence[n], energy_param[0], energy_param[1], energy_param[2]) # bp1, bp2, g_AU, g_GU, g_GC

    if is_circular:
        pred_partition = log_z.circular(g_base_pair, p.g_loop, energy_param[3], N)
        for j in range(4):
            gr = log_z.circular_derivatives(g_base_pair, p.g_loop, energy_param[3], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, g
            grad[j] += (-2)*(partition - pred_partition) * gr 
    else:
        pred_partition = log_z.linear(g_base_pair, p.g_loop, energy_param[3], N)
        for j in range(4):
            gr = log_z.linear_derivatives(g_base_pair, p.g_loop, energy_param[3], N, energy_param[j])
            grad[j] += (-2)*(partition - pred_partition) * gr
    return grad

def convergence(param, prev_param): 
    if np.linalg.norm(param - prev_param) < 1e-8:
        return True
    else:
        return False


#actual_param = [5.69, 6., 4.09, 1., -7.09]
#training_data, sequences = rsg.training_generator(50, actual_param, 10, 21)

param = np.array([4.,2.,3.,5.])
prev_param = np.zeros(4)
index = range(len(p.training_data))

m = np.zeros(4)
v = np.zeros(4)
t = 0

k = 1

K = []
iteration_param = [] # list of updated parameters after each iteration

while not convergence(param, prev_param) and k < 500:
    random.shuffle(index) # randomly shuffling indexes for each pass through the data
    print k
    for i in index:
        prev_param = np.array(param)
        t += 1
        grad = gradient(param, p.training_sequences[i], p.training_data[i])
        m = p.beta1 * m + (1 - p.beta1) * grad
        v = p.beta2 * v + (1 - p.beta2) * grad**2
        mhat = m / (1 - p.beta1**t)
        vhat = v / (1 - p.beta2**t)
        param += (-p.alpha) * mhat / (np.sqrt(vhat) + 1e-8)
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