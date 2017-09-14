# last modified: September 14, 2017
# toy model to calculate log(partition) for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# NOTE IGNORING MULTILOOPS
import time

start = time.time()

import numpy as np
import log_partition as log_z
#import random_sequence_generator as rsg
import train_toy as train
import scipy.optimize as so
import parameters as p

'''all free energy parameters in kcal/mol'''

def predicted_logpartition(energy_param, g_loop, test_sequences):
    predicted = np.zeros(len(test_sequences))
    for i in range(len(test_sequences)):
        sequence = test_sequences[i]
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
            predicted[i] = log_z.circular(g_base_pair, g_loop, energy_param[3], N) # g_base_pair, g_loop, g_stack, N
        else:
            predicted[i] = log_z.linear(g_base_pair, g_loop, energy_param[3], N) # g_base_pair, g_loop, g_stack, N
    return predicted

def frac_resid(energy_param, g_loop, sequences):
    return abs(p.test_data - predicted_logpartition(energy_param, g_loop, sequences))/p.test_data

param = np.array([4., 2., 3., 5.])

#actual_param = np.array([ 4.80680612,  1.60477761,  1.05042392, -2.38723362])
#print so.check_grad(train.cost, train.cost_gradient, actual_param, g_loop, training_data, training_sequences)
#print train.cost_gradient(actual_param, g_loop, training_data, training_sequences)
#print so.approx_fprime(actual_param, train.cost, 1e-8, g_loop, training_data, training_sequences)

optimization = so.minimize(train.cost, param, args=(p.g_loop, p.training_data, p.training_sequences), method='BFGS', jac=train.cost_gradient, tol=1e-8, callback=train.check)
#optimization = so.minimize(train.cost, param, args=(p.g_loop, p.training_data, p.training_sequences), method='BFGS', tol=1e-8, callback=train.check)

final = optimization.x
print "minimize completed"
print optimization

print '----------------------------------------\n'
print 'Test:   '+str(p.test_data)
print 'Pred:   '+str(predicted_logpartition(final, p.g_loop, p.test_sequences))
print 'fr err: '+str(frac_resid(final, p.g_loop, p.test_sequences))

end = time.time()

print '\n'+str(end-start)
