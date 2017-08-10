# last modified: August 9, 2017
# toy model to calculate partition function for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# TODO IGNORING MULTILOOPS
import numpy as np
import scipy.optimize as so
import partition as z
import random_sequence_generator as r
'''all free energy parameters in kcal/mol'''

'''
# generating training data
f = open('sequences.txt','r')

training_data = np.array(f.readline())

sequence = f.readline()[:-1]
N = len(sequence)
while sequence:
    datum = []
    datum.append(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring the '-' part of sequence
    
    # initializing base-pair free energy matrix
    # IGNORING STERIC EFFECTS
    g_base_pair = np.zeros((N,N))
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = z.free_energy_pair(sequence[m],sequence[n],5.69,6.0,4.09) # bp1, bp2, g_AU, g_GU, g_GC
    
    if is_circular:
        datum.append(z.circular_sequence_partition(g_base_pair, 1.0, -7.09, N)) # g_base_pair, g_loop, g_stack, N
    else:
        datum.append(z.linear_sequence_partition(g_base_pair, 1.0, -7.09, N)) # g_base_pair, g_loop, g_stack, N
    
    training_data.append(datum)
    sequence = f.readline()[:-1]
    N = len(sequence)

f.close()
'''

training_data = r.training_data
# sum of square residuals
# want to minimize to get least squares
# energy_param = [g_AU,g_GU,g_GC,g_loop_g_stack]
def ssr(energy_param):
    sr = 0.0
    for i in range(len(training_data)):
        sequence = training_data[i][0]
        N = len(sequence)
        x = 0.0
        
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
            x = z.circular_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
        else:
            x = z.linear_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
        
        sr += (training_data[i][1] - x)**2
    return sr

# gradient of ssr wrt energy parameters
def ssr_prime(energy_param):
    sr_p = np.zeros(5)
    for i in range(len(training_data)):
        sequence = training_data[i][0]
        N = len(sequence)
        x = 0.0
        dx = np.zeros(5)
        
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
            x = z.circular_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N 
            for j in range(5):
                dx[j] = z.circular_sequence_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, dg
        else:
            x = z.linear_sequence_partition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
            for j in range(5):
                dx[j] = z.linear_sequence_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
        sr_p += -2*(training_data[i][1] - x)*dx
    return sr_p


param = np.array([1,2,3,4,5])

#print so.fmin_bfgs(ssr,param,gtol=1e-6)
print so.fmin_bfgs(ssr,param,fprime=ssr_prime,gtol=1e-6)
