import numpy as np
import log_partition as logz
#import training_sequence_generator as tr 

# TODO need to change how training data is stored (hopefully using tables)
# returns gradient of residual squared for particular point
def gradient(energy_param,training_datum):
    grad = np.zeros(5)
    sequence = training_datum[0]
    pred_logpartition = training_datum[1]
    
    N = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring '-'

    # ignoring steric effects, generating base-pair free energy matrix
    g_base_pair = np.zeros((N,N))
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = logz.free_energy_pair(sequence[m], sequence[n], energy_param[0], energy_param[1], energy_param[2]) # bp1, bp2, g_AU, g_GU, g_GC

    if is_circular:
        for j in range(5):
            grad[j] = logz.circular_sequence_logderivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, g
        grad = grad * (-2*(pred_logpartition - logz.circular_sequence_logpartition(g_base_pair, energy_param[3], energy_param[4], N)))
    else:
        for j in range(5):
            grad[j] = logz.linear_sequence_logderivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
        grad = grad * (-2*(pred_logpartition - logz.linear_sequence_logpartition(g_base_pair, energy_param[3], energy_param[4], N)))
    
    return grad

def update_param(energy_param, learning_rate, training_datum):
    return energy_param - learning_rate * gradient(energy_param, training_datum)

def convergence(param, prev_param):
    if np.linalg.norm(param - prev_param) < 1e-8:
        return True
    else:
        return False


