# modified: November 10, 2017
# using Adam method
import numpy as np
import partition as z
#import log_partition as log_z
import base_pair_probabilities as bpp
import parameters as p
import random
import filtering as fil
import sys

import cProfile 

g_loop = p.g_loop
# returns gradient of residual squared for particular point
# NOTE new cost function: sum( (Qi - Q(xi,B))^2 ) + alpha*B^2
def gradient(param, sequence, basepairprob, bp):
    deriv_sum_sq = np.zeros(4)
    M = len(sequence)
    
    bpp_matrix = bpp.mccaskill_linear(param, sequence, M)
    residual = basepairprob - bpp_matrix[bp[0], bp[1]]
    deriv_bpp_matrix = bpp.mccaskill_linear_gradient(param, sequence, M)
    deriv_sum_sq = (-2) * residual * deriv_bpp_matrix[bp[0], bp[1]]
    prior = p.priors - param
    return deriv_sum_sq - 2 * p.w * prior #NOTE second term: prior term to deal with vanishing gradient

# NOTE will need to fix the function calls if wanting to use this section XXX
#def gradient(energy_param, sequence, partition):
#    N = len(sequence)
#
#    is_circular = False
#    if sequence[-1] == '-':
#        is_circular = True
#        N = N - 1 # ignoring '-'
#
#    prior = p.priors - energy_param
#    grad = -2 * p.w * prior
#    
#    # ignoring steric effects, generating base-pair free energy matrix
#    g_base_pair = gm.generator(sequence, energy_param[0], energy_param[1], energy_param[2], N)
#
#    if is_circular:
#        pred_partition = log_z.circular(g_base_pair, p.g_loop, energy_param[3], N)
#        for j in range(4):
#            gr = log_z.circular_derivatives(g_base_pair, p.g_loop, energy_param[3], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, g
#            grad[j] += (-2)*(partition - pred_partition) * gr 
#    else:
#        pred_partition = log_z.linear(g_base_pair, p.g_loop, energy_param[3], N)
#        for j in range(4):
#            gr = log_z.linear_derivatives(g_base_pair, p.g_loop, energy_param[3], N, energy_param[j])
#            grad[j] += (-2)*(partition - pred_partition) * gr
#    return grad

def convergence(param, prev_param): 
    if np.linalg.norm(param - prev_param) < 1e-6:
        return True
    else:
        return False

def main():

    #actual_param = [5.69, 6., 4.09, 1., -7.09]
    #training_data, sequences = rsg.training_generator(50, actual_param, 10, 21)
    num_training_examples = 50
    training_data = [] # storing bpp for closing stem of hairpin for each input sequence
    g_loop = p.g_loop
    sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
    for i in range(num_training_examples):
        training_data.append((1e-9) / fil.KDnoligand[i])
    print 'finished gathering training data'

    param = np.array([1.,2.,3.,4.])
    prev_param = np.zeros(4)
    index = range(len(training_data)) #p.training_data))

    m = np.zeros(4)
    v = np.zeros(4)
    t = 0

    k = 1

    K = []
    iteration_param = [] # list of updated parameters after each iteration

    while not convergence(param, prev_param) and k < 2000:
        random.shuffle(index) # randomly shuffling indexes for each pass through the data
        if k % 100 == 0:
            print "Gradient " + str(grad)
            print "Updated: " + str(param)
        print k
        for i in index:
            prev_param = np.array(param)
            t += 1
            grad = gradient(param, sequences[i], training_data[i], fil.closing_bp_indices[i]) #p.training_sequences[i], p.training_data[i])
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
    print iteration_param
    print param
    

profile = cProfile.Profile()
profile.enable()
main()
profile.disable()
profile.print_stats()
profile.dump_stats("%s.statout" % sys.argv[1])

