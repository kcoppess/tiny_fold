# last modified: August 10, 2017
# toy model to calculate partition function for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# TODO IGNORING MULTILOOPS
import numpy as np
import scipy.optimize as so
import partition as z
import random_sequence_generator as r

'''all free energy parameters in kcal/mol'''

def predicted_partition(energy_param):
    predicted = np.zeros(len(r.training_data))
    for i in range(len(r.training_data)):
        sequence = r.sequences[i]
        N = len(sequence)
        
        is_circular = False
        if sequence[-1] == '-':
            is_circular = True
            N = N - 1 # ignoring the '-' part of sequence
        
        # initializing base-pair free energy matrix
        # IGNORING STERIC EFFECTS
        g_base_pair = np.zeros((N,N))
        for m in range(N):
            for n in range(N):
                g_base_pair[m,n] = z.free_energy_pair(sequence[m],sequence[n], energy_param[0], energy_param[1], energy_param[2]) # bp1, bp2, g_AU, g_GU, g_GC
        
        if is_circular:
            predicted[i] = z.circular_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
        else:
            predicted[i] = z.linear_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
    return predicted

def predicted_gradient(energy_param):
    predicted_grad = np.zeros((len(r.training_data), 5))
    for i in range(len(r.training_data)):
        sequence = r.sequences[i]
        N = len(sequence)
        
        is_circular = False
        if sequence[-1] == '-':
            is_circular = True
            N = N - 1 # ignoring the '-' part of sequence
        
        # initializing base-pair free energy matrix
        # IGNORING STERIC EFFECTS
        g_base_pair = np.zeros((N,N))
        for m in range(N):
            for n in range(N):
                g_base_pair[m,n] = z.free_energy_pair(sequence[m],sequence[n], energy_param[0], energy_param[1], energy_param[2]) # bp1, bp2, g_AU, g_GU, g_GC
        
        if is_circular:
            for j in range(5):
                predicted_grad[i][j] = z.circular_sequence_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, dg
        else:
            for j in range(5):
                predicted_grad[i][j] = z.linear_sequence_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
    return predicted_grad

# sum of square residuals
# want to minimize to get least squares
# energy_param = [g_AU,g_GU,g_GC,g_loop_g_stack]
def ssr(energy_param):
    residuals = r.training_data - predicted_partition(energy_param)
    return np.dot(residuals, residuals) # TODO overflow error keeps popping up

# gradient of ssr wrt energy parameters
def ssr_prime(energy_param):
    sr_p = np.zeros(5)
    residuals = r.training_data - predicted_partition(energy_param)
    return (-2)*np.dot(residuals, predicted_gradient(energy_param))


param = np.array([2.,3.,1.,4.,-2.])

print so.minimize(ssr,param,method='Nelder-Mead',tol=1e-8)
#print so.fmin_bfgs(ssr,param,fprime=ssr_prime,gtol=1e-6)
