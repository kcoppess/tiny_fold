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

def predicted_gradient(energy_param, g_loop, sequences):
    predicted_grad = np.zeros((len(sequences), 4))
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
            for j in range(4):
                predicted_grad[i][j] = log_z.circular_derivatives(g_base_pair, g_loop, energy_param[3], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, dg
        else:
            for j in range(4):
                predicted_grad[i][j] = log_z.linear_derivatives(g_base_pair, g_loop, energy_param[3], N, energy_param[j])
    return predicted_grad

# cost function
# energy_param = [g_AU,g_GU,g_GC,g_stack]
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost(energy_param, g_loop, training_data, sequences):
    #print energy_param
    residuals = training_data - predicted(energy_param, g_loop, sequences)
    return np.dot(residuals, residuals) + p.w * np.dot(energy_param, energy_param) #NOTE second term: prior term to deal with vanishing gradient

# gradient of cost wrt energy parameters
# NOTE: g_loop is being fixed (only optimizing for 4 parameters)
def cost_gradient(energy_param, g_loop, training_data, sequences):
    residuals = training_data - predicted(energy_param, g_loop, sequences)
    grad = (-2)*np.dot(residuals, predicted_gradient(energy_param, g_loop, sequences)) + 2 * p.w * energy_param
    return grad

def check(x):
    print str(x) #+'      '+str(cost_gradient(x))
'''
#NOTE will not converge to correct parameters if one of the initial guesses is 0.0
param = np.array([1.,2.,3.,4.,5.])

#final = so.minimize(cost,param,method='BFGS',jac=cost_gradient,tol=1e-8,callback=check)
#print final
final_res = so.minimize(cost,param,method='Nelder-Mead',tol=1e-8,options={'maxiter':10000,'maxfev':10000})
final = final_res.x
print str(final_res)+'\n'
'''
