# bpp training using BFGS
# NOTE ignoring circular sequences for the moment
import numpy as np
import base_pair_probabilities as bpp
import parameters as p
import scipy.optimize as so
import time
import filtering as fil

'''all free energy parameters in kcal/mol'''

'''BPP TRAINING DATA'''
num_training_examples = 10
training_data = [] # storing bpp for closing stem of hairpin for each input sequence
#energy_param = p.energies
g_loop = p.g_loop
sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
for i in range(num_training_examples):
    training_data.append((1e-9) / fil.KDnoligand[i])
print 'finished gathering training data'

'''COST FUNCTION'''

# cost function
# energy_param = [g_AU,g_GU,g_GC,g_stack]
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost(param):
    sum_sq = 0.
    for m in range(num_training_examples):
        M = len(sequences[m])

        bp = fil.closing_bp_indices[m]
        bpp_matrix = bpp.mccaskill_linear(param, sequences[m], M)
        sq_residual = (training_data[m] - bpp_matrix[bp[0], bp[1]])**2
        sum_sq += sq_residual
    prior = p.priors - param
    return sum_sq + p.w * np.dot(prior, prior)
'''
    sum_sq = 0.
    for m in range(len(sequences)):
        bppm = training_data[m]
        M = len(sequences[m])
        if sequences[m][-1] == '-':
            M = M - 1 # ignoring the '-' part of sequence
        g_base_pair = gm.generator(sequences[m], param[0], param[1], param[2], M)
        
        residual = bppm - bpp.mccaskill_linear(param, sequences[m], M)
        square_residual = residual**2
        sum_sq += square_residual.sum()
    prior = p.priors - param
    return sum_sq + p.w * np.dot(prior, prior) #NOTE second term: prior term to deal with vanishing gradient
'''

# gradient of cost wrt energy parameters
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost_gradient(param):
    deriv_sum_sq = np.zeros(4)
    for m in range(num_training_examples):
        M = len(sequences[m])
        
        bp = fil.closing_bp_indices[m]
        bpp_matrix = bpp.mccaskill_linear(param, sequences[m], M)
        residual = training_data[m] - bpp_matrix[bp[0], bp[1]]
        deriv_bpp_matrix = bpp.mccaskill_linear_gradient(param, sequences[m], M, param[k])
        deriv_sum_sq += (-2) * residual * deriv_bpp_matrix[bp[0], bp[1]]
    prior = p.priors - param
    return deriv_sum_sq - 2 * p.w * prior #NOTE second term: prior term to deal with vanishing gradient
'''
    deriv_sum_sq = np.zeros(4)
    for m in range(len(sequences)):
        bppm = training_data[m]
        M = len(sequences[m])

        if sequences[m][-1] == '-':
            M = M - 1 # ignoring the '-' part of sequence
        
        residual = bppm - bpp.mccaskill_linear(param, sequences[m], M)
        square_residual_deriv = (-2) * residual * bpp.mccaskill_linear_gradient(param, sequences[m], M, param[k])
        deriv_sum_sq += square_residual_deriv.sum()
    prior = p.priors - param
    return deriv_sum_sq - 2 * p.w * prior #NOTE second term: prior term to deal with vanishing gradient
'''

param_iterations = []

def check(x):
    param_iterations.append(x)
    print str(cost_gradient(x))
    return

par = p.energies


start = time.clock()
#optimization = so.minimize(cost, par, method='Nelder-Mead', tol=1e-8, callback=check)
optimization = so.minimize(cost, par, method='BFGS', jac=cost_gradient, tol=1e-8, callback=check)
end = time.clock()
print 'Time: '+str(end - start)

final = optimization.x
#print "minimize completed"
print optimization

print ' ' 
print param_iterations

