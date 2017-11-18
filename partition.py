# last modified: November 10, 2017
# 
# toy model to calculate partition function for a given sequence
# based off of NUPACK pseudo-code given in Dirks & Pierce 2003
# NOTE IGNORING MULTILOOPS
#
# Q_hairpin(g_BP)
# Q_interior(g_BP, g_stack, loop_type)
#
# functions with N^3 algorithm:
#  linear(param, sequence, N)
#  linear_gradient(param, sequence, N, g)
#  linear_derivatives(param, sequence, N, g)
#  linear_derivatives_over_val(param, sequence, N, g)
#  circular(param, sequence, N)
#  circular_gradient(param, sequence, N)
#  circular_derivatives(param, sequence, N, g)
#  circular_derivatives_over_val(param, sequence, N, g)
#
# analytical derivatives
#
import numpy as np
import parameters as p
import g_matrix as gm

'''all free energy parameters in kcal/mol'''

R = p.R # kcal/K/mol universal gas constant
T = p.T # K temperature (standard state - room temp)
invRT = 1.0 / (R*T)
g_loop = p.g_loop

# partition contribution from hairpin closed with base pair BP
def Q_hairpin(g_BP):
    return np.exp(-(g_BP + g_loop)/(R*T))

# partition contribution from interior loop closed with base pair BP
# avoids double counting free energy from base pair formation
def Q_interior(g_BP, g_stack, loop_type):
    if loop_type == 's':  # if stacked base pair
        return np.exp(-(g_BP + g_stack)*invRT)
    elif loop_type == 'l':  # if loop
        return np.exp(-(g_BP + g_loop)*invRT)
    else:
        print str(loop_type)+": Invalid type of interior loop"
        return 0.0

def linear(param, sequence, N):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    Qs = np.zeros((N,N))    
    # calculation of partition function
    for ll in range(1,N+1): # iterating over all subsequence lengths
        for ii in range(0,N-ll+1): # iterating over all starting positions for subsequences
            jj = ii + ll - 1 # ending position for subsequence
            # Qb recursion
            if jj-ii > 3 and g_base_pair[ii,jj]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                Qb[ii,jj] = Q_hairpin(g_base_pair[ii,jj])
            else: # no hairpin possible
                Qb[ii,jj] = 0.0
            for dd in range(ii+1,jj-4): # iterate over all possible rightmost pairs
                for ee in range(dd+4,jj): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if g_base_pair[ii,jj] and g_base_pair[dd,ee]: # possible for both base pairs to form
                        if ii+1 == dd and ee+1 == jj: # if stacked
                            interior_loop_type = 's'
                        else: # if loop
                            interior_loop_type = 'l'
                        Qb[ii,jj] += Qb[dd,ee] * Q_interior(g_base_pair[ii,jj], g_stack, interior_loop_type)
                    else: # no interior loop possible (one or both base pairs can't form)
			pass
            # Qs recursion
            for dd in range(ii+4, jj+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[ii,jj] += Qb[ii,dd]
            # Q recursion
            Q[ii,jj] = 1.0
            for dd in range(ii,jj-3): # iterating over all possible rightmost pairs
                if dd == 0: # to deal with issue of wrapping around in the last iteration
                    Q[ii,jj] += Qs[dd,jj]
                else:
                    Q[ii,jj] += Q[ii,dd-1]*Qs[dd,jj]
    return Q[0,N-1]


def linear_derivatives(param, sequence, N, g): # g: energy parameter differentiating w.r.t.
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)

    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    dQ = np.zeros((N,N)) # stores derivatives

    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N)) # stores bound derivatives

    Qs = np.zeros((N,N))
    dQs = np.zeros((N,N))

    # LINEAR SEQUENCE calculation of partition function and its derivative
    for ll in range(1,N+1): # iterating over all subsequence lengths
        for ii in range(0,N-ll+1): # iterating over all starting positions for subsequences
            jj = ii + ll - 1 # ending position for subsequence
            # Qb recursion
            if jj-ii > 3 and g_base_pair[ii,jj]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                q_hairpin_ij = Q_hairpin(g_base_pair[ii,jj])
                Qb[ii,jj] = q_hairpin_ij
                if g == g_base_pair[ii,jj]:  # differentiating wrt current base pair parameter
                    dQb[ii,jj] = -q_hairpin_ij * invRT #(Q_hairpin(g_base_pair[i,j] + h, g_loop) - q_hairpin_ij)/h
            else: # no hairpin possible
                Qb[ii,jj] = 0.0
                dQb[ii,jj] = 0.0
            for dd in range(ii+1,jj-4): # iterate over all possible rightmost pairs
                for ee in range(dd+4,jj): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = '' #g_interior = 0.0
                    if g_base_pair[ii,jj] and g_base_pair[dd,ee]: # possible for both base pairs to form
                        if ii+1 == dd and ee+1 == jj: # if stacked
                            interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack # avoids double counting free energy from base pair formation
                        else: # if loop
                            interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                        q_int_ij = Q_interior(g_base_pair[ii,jj], g_stack, interior_loop_type)
                        Qb[ii,jj] += Qb[dd,ee] * q_int_ij
                        dQ_int = 0.0 # derivative of Q_interior
                        if g == g_base_pair[ii,jj]:
                            dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - q_int_ij)/h
                        elif g == g_stack and interior_loop_type == 's':
                            dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - q_int_ij)/h
                        dQb[ii,jj] += dQ_int * Qb[dd,ee] + q_int_ij * dQb[dd,ee]
                    else: # no interior loop possible (one or both base pairs can't form)
                        pass
            # Qs recursion
            for dd in range(ii+4, jj+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[ii,jj] += Qb[ii,dd]
                dQs[ii,jj] += dQb[ii,dd]
            # Q recursion
            Q[ii,jj] = 1.0
            for dd in range(ii,jj-3): # iterating over all possible rightmost pairs
                if dd == 0: # to deal with issue of wrapping around in the last iteration
                    Q[ii,jj] += Qs[dd,jj]
                    dQ[ii,jj] += dQs[dd,jj]
                else:
                    Q[ii,jj] += Q[ii,dd-1]*Qs[dd,jj]
                    dQ[ii,jj] += dQ[ii,dd-1]*Qs[dd,jj] + Q[ii,dd-1]*dQs[dd,jj]
    return dQ[0,N-1]

def linear_gradient(param, sequence, N):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    dQ = np.zeros((N,N,4)) # stores gradients

    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N,4)) # stores bound derivatives

    Qs = np.zeros((N,N))
    dQs = np.zeros((N,N,4))

    # LINEAR SEQUENCE calculation of partition function and its derivative
    for ll in range(1,N+1): # iterating over all subsequence lengths
        for ii in range(0,N-ll+1): # iterating over all starting positions for subsequences
            jj = ii + ll - 1 # ending position for subsequence
            # Qb recursion
            energy_index, = np.where(param == g_base_pair[ii,jj])
            if jj-ii > 3 and g_base_pair[ii,jj]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                q_hairpin_ij = Q_hairpin(g_base_pair[ii,jj])
                Qb[ii,jj] = q_hairpin_ij
                dQb[ii,jj, energy_index] = -q_hairpin_ij * invRT #(Q_hairpin(g_base_pair[i,j] + h, g_loop) - q_hairpin_ij)/h
            else: # no hairpin possible
                pass
            for dd in range(ii+1,jj-4): # iterate over all possible rightmost pairs
                for ee in range(dd+4,jj): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = '' #g_interior = 0.0
                    if g_base_pair[ii,jj] and g_base_pair[dd,ee]: # possible for both base pairs to form
                        if ii+1 == dd and ee+1 == jj: # if stacked
                            interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack # avoids double counting free energy from base pair formation
                        else: # if loop
                            interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                        q_int_ij = Q_interior(g_base_pair[ii,jj], g_stack, interior_loop_type)
                        Qb[ii,jj] += Qb[dd,ee] * q_int_ij
                        dQ_int = np.zeros(4) # derivative of Q_interior
                        dQ_int[energy_index] = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - q_int_ij)/h
                        if interior_loop_type == 's':
                            dQ_int[3] = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - q_int_ij)/h
                        dQb[ii,jj] += dQ_int * Qb[dd,ee] + q_int_ij * dQb[dd,ee]
                    else: # no interior loop possible (one or both base pairs can't form)
                        pass
            # Qs recursion
            for dd in range(ii+4, jj+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[ii,jj] += Qb[ii,dd]
                dQs[ii,jj] += dQb[ii,dd]
            # Q recursion
            Q[ii,jj] = 1.0
            for dd in range(ii,jj-3): # iterating over all possible rightmost pairs
                if dd == 0: # to deal with issue of wrapping around in the last iteration
                    Q[ii,jj] += Qs[dd,jj]
                    dQ[ii,jj] += dQs[dd,jj]
                else:
                    Q[ii,jj] += Q[ii,dd-1]*Qs[dd,jj]
                    dQ[ii,jj] += dQ[ii,dd-1]*Qs[dd,jj] + Q[ii,dd-1]*dQs[dd,jj]
    return dQ[0,N-1]

def linear_derivatives_over_val(param, sequence, N, g):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    dQ = np.zeros((N,N)) # stores derivatives

    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N)) # stores bound derivatives

    Qs = np.zeros((N,N))
    dQs = np.zeros((N,N))

    # LINEAR SEQUENCE calculation of partition function and its derivative
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                q_hairpin_ij = Q_hairpin(g_base_pair[i,j])
                Qb[i,j] = q_hairpin_ij
                if g == g_base_pair[i,j]:  # differentiating wrt current base pair parameter
                    dQb[i,j] = -q_hairpin_ij * invRT #(Q_hairpin(g_base_pair[i,j] + h, g_loop) - q_hairpin_ij)/h
            else: # no hairpin possible
                Qb[i,j] = 0.0
                dQb[i,j] = 0.0
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = '' #g_interior = 0.0
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack # avoids double counting free energy from base pair formation
                        else: # if loop
                            interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                        q_int_ij = Q_interior(g_base_pair[i,j], g_stack, interior_loop_type)
                        Qb[i,j] += Qb[d,e] * q_int_ij
                        dQ_int = 0.0 # derivative of Q_interior
                        if g == g_base_pair[i,j]:
                            dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - q_int_ij)/h
                        elif g == g_stack and interior_loop_type == 's':
                            dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - q_int_ij)/h
                        dQb[i,j] += dQ_int * Qb[d,e] + q_int_ij * dQb[d,e]
                    else: # no interior loop possible (one or both base pairs can't form)
                        pass
            
            # Qs recursion
            for d in range(i+4, j+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[i,j] += Qb[i,d]
                dQs[i,j] += dQb[i,d]
            # Q recursion
            Q[i,j] = 1.0
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                if d == 0: # to deal with issue of wrapping around in the last iteration
                    Q[i,j] += Qs[d,j]
                    dQ[i,j] += dQs[d,j]
                else:
                    Q[i,j] += Q[i,d-1]*Qs[d,j]
                    dQ[i,j] += dQ[i,d-1]*Qs[d,j] + Q[i,d-1]*dQs[d,j]
    return dQ[0,N-1]/Q[0,N-1]

def circular(param, sequence, N):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0

    Qb = np.zeros((N,N))
    Qs = np.zeros((N,N))

    exp_neg_gloop_over_RT = np.exp(-g_loop*invRT)
    exp_neg_gstack_gloop_over_RT = np.exp(-(g_stack-g_loop)*invRT)

    for l in range(1, N+1): #length of subsequence
        for i in range(0, N-l+1): # start of subsequence
            j = i + l - 1
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                Qb[i,j] = Q_hairpin(g_base_pair[i,j])
            else:  # hairpin not possible
                Qb[i,j] = 0.0
            
            for d in range(i+1,j-4):
                for e in range(d+4,j):
                    interior_loop_type = ''
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]: # interior loop possible
                            if i+1 == d and e+1 == j: #stacked
                                interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack
                            else: #loop
                                interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                            Qb[i,j] += Qb[d,e] * Q_interior(g_base_pair[i,j], g_stack, interior_loop_type)
                        else: # interior loop not possible
                            pass
                    else:
                        pass
            # Qs recursion
            for d in range(i+4, j+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[i,j] += Qb[i,d]
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    if d == 0:
                        Q[i,j] += Qs[0,j]*exp_neg_gloop_over_RT
                    else:
                        if Qb[0,d-1] and Qb[d,N-1]: # to account for stacked pair forming when chain is closed
                            Q[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                        else: # to account for interior loop forming when chain is closed
                            Q[i,j] += Q[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT
            else:
                for d in range(i,j-3):
                    if d == 0:
                        Q[i,j] += Qs[d,j]
                    else:
                        Q[i,j] += Q[i,d-1]*Qs[d,j]
    return Q[0,N-1]

def circular_gradient(param, sequence, N):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    dQ = np.zeros((N,N,4)) # stores gradients

    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N,4)) # stores bound derivatives

    Qs = np.zeros((N,N))
    dQs = np.zeros((N,N,4))

    # CIRCULAR SEQUENCE calculation of partition function and its derivative
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            energy_index, = np.where(param == g_base_pair[i,j])
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                q_hairpin_ij = Q_hairpin(g_base_pair[i,j])
                Qb[i,j] = q_hairpin_ij
                dQb[i,j, energy_index] = -q_hairpin_ij * invRT #(Q_hairpin(g_base_pair[i,j] + h, g_loop) - q_hairpin_ij)/h
            else: # no hairpin possible
                pass
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = '' #g_interior = 0.0
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                            if i+1 == d and e+1 == j: # if stacked
                                interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack # avoids double counting free energy from base pair formation
                            else: # if loop
                                interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                            q_int_ij = Q_interior(g_base_pair[i,j], g_stack, interior_loop_type)
                            Qb[i,j] += Qb[d,e] * q_int_ij
                            dQ_int = np.zeros(4) # derivative of Q_interior
                            dQ_int[energy_index] = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - q_int_ij)/h
                            if interior_loop_type == 's':
                                dQ_int[3] = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - q_int_ij)/h
                            dQb[i,j] += dQ_int * Qb[d,e] + q_int_ij * dQb[d,e]
                        else: # no interior loop possible (one or both base pairs can't form)
                            pass
                    else:
                        pass
            # Qs recursion
            for d in range(i+4, j+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[i,j] += Qb[i,d]
                dQs[i,j] += dQb[i,d]
            # Q recursion
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    exp_neg_gloop_over_RT = np.exp(-g_loop*invRT)
                    if d == 0:
                        Q[i,j] += Qs[0,j]* exp_neg_gloop_over_RT
                        dQ[i,j] += dQs[0,j] * exp_neg_gloop_over_RT
                    else:
                        if Qb[0,d-1] and Qb[d,N-1]: # to account for stacked pair forming when chain is closed
                            Q[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                            dQ[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*dQs[d,N-1]*exp_neg_gloop_over_RT
                            dQ[i,j] += (dQ[0,d-1] + dQs[0,d-1]*(exp_neg_gstack_gloop_over_RT -1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                            dQ[i,j,3] += Qs[0,d-1]*(-exp_neg_gstack_gloop_over_RT)*Qs[d,N-1]*exp_neg_gloop_over_RT
                        else: # to account for interior loop forming when chain is closed
                            Q[i,j] += Q[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT
                            dQ[i,j] += dQ[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT + Q[i,d-1]*dQs[d,j]*exp_neg_gloop_over_RT
            else:
                for d in range(i,j-3):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i,j] += Qs[d,j]
                        dQ[i,j] += dQs[d,j]
                    else:
                        Q[i,j] += Q[i,d-1]*Qs[d,j]
                        dQ[i,j] += dQ[i,d-1]*Qs[d,j] + Q[i,d-1]*dQs[d,j]
    return dQ[0,N-1]

def circular_derivatives(param, sequence, N, g):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)

    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0

    dQ = np.zeros((N,N))

    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N))

    Qs = np.zeros((N,N))
    dQs = np.zeros((N,N))
       
    exp_neg_gstack_gloop_over_RT = np.exp(-(g_stack-g_loop)*invRT)

    for l in range(1, N+1): #length of subsequence
        for i in range(0, N-l+1): # start of subsequence
            j = i + l - 1
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                q_hairpin_ij = Q_hairpin(g_base_pair[i,j])
                Qb[i,j] = q_hairpin_ij
                if g == g_base_pair[i,j]:
                    dQb[i,j] = -q_hairpin_ij * invRT #(Q_hairpin(g_base_pair[i,j] + h, g_loop) - q_hairpin_ij)/h
            else:  # hairpin not possible
				Qb[i,j] = 0.0
            for d in range(i+1,j-4):
                for e in range(d+4,j):
                    interior_loop_type = ''
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]: # interior loop possible
                            if i+1 == d and e+1 == j: #stacked
                                interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack
                            else: #loop
                                interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                            q_int_ij = Q_interior(g_base_pair[i,j], g_stack, interior_loop_type)
                            Qb[i,j] += Qb[d,e] * q_int_ij
                            dQ_int = 0.0
                            if g == g_base_pair[i,j]:
                                dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - q_int_ij)/h
                            elif g == g_stack and interior_loop_type == 's':
                                dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - q_int_ij)/h
                            dQb[i,j] += dQ_int * Qb[d,e] + q_int_ij * dQb[d,e]
                        else: # interior loop not possible
                            pass
                    else:
                        pass
            for d in range(i+4, j+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[i,j] += Qb[i,d]
                dQs[i,j] += dQb[i,d]
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    exp_neg_gloop_over_RT = np.exp(-g_loop*invRT)
                    if d == 0:
                        Q[i,j] += Qs[0,j]* exp_neg_gloop_over_RT
                        dQ[i,j] += dQs[0,j] * exp_neg_gloop_over_RT
                    else:
                        if Qb[0,d-1] and Qb[d,N-1]: # to account for stacked pair forming when chain is closed
                            Q[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                            if g == g_stack:
                                dQ[i,j] += (dQ[0,d-1] + dQs[0,d-1]*(exp_neg_gstack_gloop_over_RT -1) + Qs[0,d-1]*(-exp_neg_gstack_gloop_over_RT))*Qs[d,N-1]*exp_neg_gloop_over_RT
                                dQ[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*dQs[d,N-1]*exp_neg_gloop_over_RT
                            else:
                                dQ[i,j] += (dQ[0,d-1] + dQs[0,d-1]*(exp_neg_gstack_gloop_over_RT -1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                                dQ[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*dQs[d,N-1]*exp_neg_gloop_over_RT
                        else: # to account for interior loop forming when chain is closed
                            Q[i,j] += Q[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT
                            dQ[i,j] += dQ[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT + Q[i,d-1]*dQs[d,j]*exp_neg_gloop_over_RT
            else:
                for d in range(i,j-3):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i,j] += Qs[d,j]
                        dQ[i,j] += dQs[d,j]
                    else:
                        Q[i,j] += Q[i,d-1]*Qs[d,j]
                        dQ[i,j] += dQ[i,d-1]*Qs[d,j] + Q[i,d-1]*dQs[d,j]
    return dQ[0,N-1]

def circular_derivatives_over_val(param, sequence, N, g):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)

    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    dQ = np.zeros((N,N))

    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N))

    Qs = np.zeros((N,N))
    dQs = np.zeros((N,N))

    exp_neg_gstack_gloop_over_RT = np.exp(-(g_stack-g_loop)*invRT)

    for l in range(1, N+1): #length of subsequence
        for i in range(0, N-l+1): # start of subsequence
            j = i + l - 1
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                q_hairpin_ij = Q_hairpin(g_base_pair[i,j])
                Qb[i,j] = q_hairpin_ij
                if g == g_base_pair[i,j]:
                    dQb[i,j] = -q_hairpin_ij * invRT #(Q_hairpin(g_base_pair[i,j] + h, g_loop) - q_hairpin_ij)/h
            else:  # hairpin not possible
                Qb[i,j] = 0.0
                for d in range(i+1,j-4):
                    for e in range(d+4,j):
                        interior_loop_type = ''
                        if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                            if g_base_pair[i,j] and g_base_pair[d,e]: # interior loop possible
                                if i+1 == d and e+1 == j: #stacked
                                    interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack
                                else: #loop
                                    interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                                q_int_ij = Q_interior(g_base_pair[i,j], g_stack, interior_loop_type)
                                Qb[i,j] += Qb[d,e] * q_int_ij
                                dQ_int = 0.0
                                if g == g_base_pair[i,j]:
                                    dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - q_int_ij)/h
                                elif g == g_stack and interior_loop_type == 's':
                                    dQ_int = -q_int_ij * invRT #(Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - q_int_ij)/h
                                dQb[i,j] += dQ_int * Qb[d,e] + q_int_ij * dQb[d,e]
                            else: # interior loop not possible
                                pass
                        else:
                            pass
            for d in range(i+4, j+1): # iterate over all rightmost pairs with base i (beginning of subsequence)
                Qs[i,j] += Qb[i,d]
                dQs[i,j] += dQb[i,d]
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    exp_neg_gloop_over_RT = np.exp(-g_loop*invRT)
                    if d == 0:
                        Q[i,j] += Qs[0,j]* exp_neg_gloop_over_RT
                        dQ[i,j] += dQs[0,j] * exp_neg_gloop_over_RT
                    else:
                        if Qb[0,d-1] and Qb[d,N-1]: # to account for stacked pair forming when chain is closed
                            Q[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                            if g == g_stack:
                                dQ[i,j] += (dQ[0,d-1] + dQs[0,d-1]*(exp_neg_gstack_gloop_over_RT -1) + Qs[0,d-1]*(-exp_neg_gstack_gloop_over_RT)*Qs[d,N-1]*exp_neg_gloop_over_RT)
                                dQ[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*dQs[d,N-1]*exp_neg_gloop_over_RT
                            else:
                                dQ[i,j] += (dQ[0,d-1] + dQs[0,d-1]*(exp_neg_gstack_gloop_over_RT -1))*Qs[d,N-1]*exp_neg_gloop_over_RT
                                dQ[i,j] += (Q[0,d-1] + Qs[0,d-1]*(exp_neg_gstack_gloop_over_RT - 1))*dQs[d,N-1]*exp_neg_gloop_over_RT
                        else: # to account for interior loop forming when chain is closed
                            Q[i,j] += Q[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT
                            dQ[i,j] += dQ[i,d-1]*Qs[d,j]*exp_neg_gloop_over_RT + Q[i,d-1]*dQs[d,j]*exp_neg_gloop_over_RT
            else:
                for d in range(i,j-3):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i,j] += Qs[d,j]
                        dQ[i,j] += dQs[d,j]
                    else:
                        Q[i,j] += Q[i,d-1]*Qs[d,j]
                        dQ[i,j] += dQ[i,d-1]*Qs[d,j] + Q[i,d-1]*dQs[d,j]
    return dQ[0,N-1]/Q[0,N-1]

