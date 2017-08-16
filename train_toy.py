# last modified: August 16, 2017
# training toy model to calculate log(partition) for a given sequence
import numpy as np
import scipy.optimize as so
import log_partition as log_z
import training_sequence_generator as tr

'''all free energy parameters in kcal/mol'''

def predicted_logpartition(energy_param):
    predicted = np.zeros(len(tr.training_data))
    for i in range(len(tr.training_data)):
        sequence = tr.sequences[i]
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

def predicted_gradient(energy_param):
    predicted_grad = np.zeros((len(tr.training_data), 5))
    for i in range(len(tr.training_data)):
        sequence = tr.sequences[i]
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
                    predicted_grad[i][j] = (-1/(log_z.R*log_z.T)) + log_z.circular_sequence_logderivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
                else:
                    predicted_grad[i][j] = log_z.circular_sequence_logderivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j]) # g_base_pair, g_loop, g_stack, N, dg
        else:
            for j in range(5):
                if j == 3:
                    predicted_grad[i][j] = (-1/(log_z.R*log_z.T)) + log_z.linear_sequence_logderivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
                else:
                    predicted_grad[i][j] = log_z.linear_sequence_logderivatives(g_base_pair, energy_param[3], energy_param[4], N, energy_param[j])
    return predicted_grad

# sum of square residuals
# want to minimize to get least squares
# energy_param = [g_AU,g_GU,g_GC,g_loop_g_stack]
def ssr(energy_param):
    #NOTE -g_loop/R*T put in to try to break degeneracy
    residuals = (tr.training_data - energy_param[3]/(log_z.R*log_z.T)) - (predicted_logpartition(energy_param) - tr.g_loop/(log_z.R*log_z.T))
    return np.dot(residuals, residuals) #+ (100000000000000*(tr.g_loop - energy_param[3]))**2 #NOTE second term to help deal with degeneracy

def resid(energy_param):
    return tr.training_data-predicted_logpartition(energy_param)

# gradient of ssr wrt energy parameters
def ssr_prime(energy_param):
    sr_p = np.zeros(5)
    residuals = (tr.training_data - energy_param[3]/(log_z.R*log_z.T)) - (predicted_logpartition(energy_param) - tr.g_loop/(log_z.R*log_z.T))
    return (-2)*np.dot(residuals, predicted_gradient(energy_param))

def check(x):
    print ssr(x)
#    print '{0: 3.6f}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}    {5: 3.6f}'.format(x[0], x[1], x[2], x[3], x[4], predicted_partition(x))

#NOTE will not converge to correct parameters if one of the initial guesses is 0.0
param = np.array([1.,2.,3.,4.,5.])

final = so.fmin_bfgs(ssr,param,fprime=ssr_prime,gtol=1e-8,epsilon=log_z.h)
#final_res = so.minimize(ssr,param,method='Nelder-Mead',tol=1e-8,options={'maxiter':10000,'maxfev':10000})
#final = final_res.x
print str(final)+'\n'
