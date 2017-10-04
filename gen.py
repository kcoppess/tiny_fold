# ad hoc implementation of training on bpp data
# NOTE ignoring circular sequences for the moment
import numpy as np
import base_pair_probabilities as bpp
import parameters as p
import g_matrix as gm
import scipy.optimize as so
import time

'''all free energy parameters in kcal/mol'''

'''GENERATING BPP TRAINING DATA'''
training_data = [] # storing bpp matrices for each input sequence
energy_param = p.energies
g_loop = p.g_loop
sequences = p.training_sequences[:10]
for sequence in sequences:
    N = len(sequence)

    if sequence[-1] == '-':
        N = N - 1 # ignoring the '-' part of sequence
    
    g_base_pair = gm.generator(sequence, energy_param[0], energy_param[1], energy_param[2], N)
    
    training_data.append(bpp.mccaskill_linear(g_base_pair, g_loop, energy_param[3], N))

'''COST FUNCTION'''

# cost function
# energy_param = [g_AU,g_GU,g_GC,g_stack]
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost(param):
    sum_sq = 0.
    for m in range(len(sequences)):
        bppm = training_data[m]
        M = len(bppm)
        g_base_pair = gm.generator(sequences[m], param[0], param[1], param[2], M)
        for i in range(M):
            for j in range(i+4,M):
                if g_base_pair[i,j]: # only training on nonzero entries (hopefully this will speed things up)
                    sum_sq += (bppm[i,j] - bpp.flag_linear(i, j, g_base_pair, g_loop, param[3], M))**2
    prior = p.priors - param
    return sum_sq + p.w * np.dot(prior, prior) #NOTE second term: prior term to deal with vanishing gradient


# gradient of cost wrt energy parameters
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost_gradient(param):
    sum_sq = []
    grad = np.zeros(4)
    for m in range(len(sequences)):
        bppm = training_data[m]
        M = len(bppm)
        g_base_pair = gm.generator(sequences[m], param[0], param[1], param[2], M)
        bpp_mat = bpp.mccaskill_linear(g_base_pair, g_loop, param[3], M)
        #print 'finished making bpp_mat and starting to run through g_base_pair '+str(time.clock())
        for i in range(M):
            for j in range(i+4,M):
                #print 'starting one sum '+str(time.clock())
                if g_base_pair[i,j]: # only training on nonzero entries (hopefully this will speed things up)
                    for k in range(4):
                        grad[k] = (-2)*(bppm[i,j] - bpp_mat[i,j]) * bpp.flag_linear_derivatives(i, j, g_base_pair, g_loop, param[3], M, param[k])
                    sum_sq.append(grad)
                #print 'finishing one sum '+str(time.clock())
        #print 'finished running through g_base_pair '+str(time.clock())
        #print '-------------------------------'
    prior = p.priors - param
    return np.sum(np.array(sum_sq), axis=0) - 2 * p.w * prior #NOTE second term: prior term to deal with vanishing gradient

def check(x):
    print str(cost_gradient(x))
    return

par = np.array([1., 2., 3., 4.])

#actual_param = np.array([ 4.80680612,  1.60477761,  1.05042392, -2.38723362])
#print so.check_grad(train.cost, train.cost_gradient, actual_param, g_loop, training_data, training_sequences)
#print train.cost_gradient(actual_param, g_loop, training_data, training_sequences)
#start = time.clock()
#psajkfldaj = so.approx_fprime(par, cost, 1e-8)
#print psajkfldaj
#end = time.clock()
#print end-start

#start = time.clock()
#pfdjsklaf = cost_gradient(par)
#print pfdjsklaf
#end = time.clock()
#print end-start

start = time.clock()
#optimization = so.minimize(cost, par, method='Nelder-Mead', tol=1e-8, callback=check)
optimization = so.minimize(cost, par, method='BFGS', jac=cost_gradient, tol=1e-8, callback=check)
end = time.clock()
print 'Time: '+str(end - start)

#final = optimization.x
#print "minimize completed"
print optimization


