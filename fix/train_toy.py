# last modified: August 27, 2017
# training toy model to calculate log(partition) for a given sequence
import numpy as np
import scipy.optimize as so
import log_partition as log_z
#import random_sequence_generator as rsg

alpha = 0.000001 # weight on prior term in cost

'''all free energy parameters in kcal/mol'''

def predicted(energy_param, sequences):
    predicted = np.zeros(len(sequences))
    for i in range(len(sequences)):
        sequence = sequences[i]
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
            predicted[i] = log_z.circular(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
        else:
            predicted[i] = log_z.linear(g_base_pair, energy_param[3], energy_param[4], N) # g_base_pair, g_loop, g_stack, N
    return predicted

def predicted_gradient(energy_param, sequences):
    predicted_grad = np.zeros((len(sequences), 5))
    for i in range(len(sequences)):
        sequence = sequences[i]
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
            for j in range(5):
                if j == 3:
                    predicted_grad[i][j] = (-1/(log_z.R*log_z.T)) + log_z.circular_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
                else:
                    predicted_grad[i][j] = log_z.circular_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, dg
        else:
            for j in range(5):
                if j == 3:
                    predicted_grad[i][j] = (-1/(log_z.R*log_z.T)) + log_z.linear_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
                else:
                    predicted_grad[i][j] = log_z.linear_derivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
    return predicted_grad

# cost function
# energy_param = [g_AU,g_GU,g_GC,g_loop_g_stack]
# g_loop = here actual g_loop TODO need to figure out way to break degeneracy in different way
def cost(energy_param, training_data, sequences, g_loop):
    print energy_param
    #NOTE -g_loop/R*T put in to try to break degeneracy
    residuals = (training_data - g_loop/(log_z.R*log_z.T)) - (predicted(energy_param, sequences) - energy_param[3]/(log_z.R*log_z.T))
    #residuals = tr.training_data - predicted_partition(energy_param)
    return np.dot(residuals, residuals) + alpha * np.dot(energy_param, energy_param) #NOTE second term: prior term to deal with vanishing gradient

# gradient of cost wrt energy parameters
#XXX not gradient for cost function with the degeneracy-breaking factor
def cost_gradient(energy_param, training_data, sequences, g_loop):
    residuals = (training_data - g_loop/(log_z.R*log_z.T)) - (predicted(energy_param, sequences) - energy_param[3]/(log_z.R*log_z.T))
    #residuals = tr.training_data - predicted_partition(energy_param)
    return (-2)*np.dot(residuals, predicted_gradient(energy_param, sequences) + np.array([0., 0., 0., -1./(log_z.R*log_z.T), 0.])) + 2 * alpha * energy_param

def check(x):
    print str(x)+'      '+str(cost_gradient(x))
'''
#NOTE will not converge to correct parameters if one of the initial guesses is 0.0
param = np.array([1.,2.,3.,4.,5.])

#final = so.minimize(cost,param,method='BFGS',jac=cost_gradient,tol=1e-8,callback=check)
#print final
final_res = so.minimize(cost,param,method='Nelder-Mead',tol=1e-8,options={'maxiter':10000,'maxfev':10000})
final = final_res.x
print str(final_res)+'\n'
'''