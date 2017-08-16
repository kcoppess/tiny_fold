# last modified: August 16, 2017
# toy model to calculate log(partition) for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# NOTE IGNORING MULTILOOPS
import time

start = time.time()

import numpy as np
import log_partition as log_z
import test_sequence_generator as test
import train_toy as train

'''all free energy parameters in kcal/mol'''

def predicted_logpartition(energy_param):
    predicted = np.zeros(len(test.test_data))
    for i in range(len(test.test_data)):
        sequence = test.sequences[i]
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
                g_base_pair[m,n] = log_z.free_energy_pair(sequence[m],sequence[n], energy_param[0], energy_param[1], energy_param[2]) # bp1, bp2, g_AU, g_GU, g_GC
        
        if is_circular:
            predicted[i] = log_z.circular_sequence_logpartition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
        else:
            predicted[i] = log_z.linear_sequence_logpartition(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
    return predicted

def frac_resid(energy_param):
    return abs(test.test_data - predicted_logpartition(energy_param))/test.test_data

model_param = train.final

print '----------------------------------------\n'
print 'Test:   '+str(test.test_data)
print 'Pred:   '+str(predicted_logpartition(model_param))
print 'fr err: '+str(frac_resid(model_param))

end = time.time()

print '\n'+str(end-start)
