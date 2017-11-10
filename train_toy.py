# last modified: September 14, 2017
# training toy model to calculate log(partition) for a given sequence
import numpy as np
import scipy.optimize as so
import log_partition as log_z
import parameters as p

'''all free energy parameters in kcal/mol'''

def predicted(energy_param, g_loop, sequences):
    predicted = np.zeros(len(sequences))
    for i in range(len(sequences)):
        sequence = sequences[i]
        N = len(sequence)
        
        is_circular = False
        if sequence[-1] == '-':
            is_circular = True
            N = N - 1 # ignoring the '-' part of sequence
        
        if is_circular:
            predicted[i] = log_z.circular(energy_param, sequence, N)
        else:
            predicted[i] = log_z.linear(energy_param, sequence, N)
    return predicted

def predicted_gradient(energy_param, g_loop, sequences):
    predicted_grad = np.zeros((len(sequences), 4))
    for i in range(len(sequences)):
        sequence = sequences[i]
        N = len(sequence)
        
        is_circular = False
        if sequence[-1] == '-':
            is_circular = True
            N = N - 1 # ignoring the '-' part of sequence
        
        if is_circular:
            predicted_grad[i] = log_z.circular_derivatives(energy_param, sequence, N)
        else:
            predicted_grad[i] = log_z.linear_derivatives(energy_param, sequence, N)
    return predicted_grad

# cost function
# energy_param = [g_AU,g_GU,g_GC,g_stack]
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost(energy_param, g_loop, training_data, sequences):
    #print energy_param
    residuals = training_data - predicted(energy_param, g_loop, sequences)
    prior = p.priors - energy_param
    return np.dot(residuals, residuals) + p.w * np.dot(prior, prior) #NOTE second term: prior term to deal with vanishing gradient

# gradient of cost wrt energy parameters
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost_gradient(energy_param, g_loop, training_data, sequences):
    residuals = training_data - predicted(energy_param, g_loop, sequences)
    prior = p.priors - energy_param
    grad = (-2)*np.dot(residuals, predicted_gradient(energy_param, g_loop, sequences)) - 2 * p.w * prior
    return grad

def check(x):
    print str(x) #+'      '+str(cost_gradient(x))
