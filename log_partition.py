# last modified: November 10, 2017
# toy model to calculate log(partition function) for a given sequence
# based off of NUPACK pseudo-code given in Dirks & Pierce 2003
# NOTE IGNORING MULTILOOPS
#
# g_interior(g_BP, g_stack, loop_type)
#
# functions w/ N^3 algorithm:
#  linear(param, sequence, N)
#  linear_derivatives(param, sequence, N, g)
#  circular(param, sequence, N)
#  circular_derivatives(param, sequence, N, g)
#
# analytical derivatives
# 
import numpy as np
import scipy.misc as sm
import partition as z
import parameters as p
import g_matrix as gm
'''all free energy parameters in kcal/mol'''

R = p.R # kcal/K/mol universal gas constant
T = p.T # K temperature (standard state - room temp)
invRT = 1.0/(R*T)
g_loop = p.g_loop

# free energy contribution from interior loop closed with base pair BP
# avoids double counting free energy from base pair formation
def g_interior(g_BP, g_stack, loop_type):
    if loop_type == 's':  # if stacked base pair
        return g_BP + g_stack
    elif loop_type == 'l':  # if loop
        return g_BP + g_loop
    else:
        print str(loop_type)+": Invalid type of interior loop"
        return 0.0

def linear(param, sequence, N):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)

    # initializing general partition matrix
    Q = [[[0] for _ in range(N)] for _ in range(N)] # stores a,b,c, etc in exp(a) + exp(b) + ...
    
    # initializing bound partition matrix
    Qb = [[[] for _ in range(N)] for _ in range(N)]
    Qs = [[[] for _ in range(N)] for _ in range(N)]

    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                Qb[i][j].append(-(g_base_pair[i,j] + g_loop)*invRT)
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            interior_loop_type = 's'
                        else: # if loop
                            interior_loop_type = 'l'
                        Qb[i][j] += [ qb - (g_interior(g_base_pair[i,j], g_stack, interior_loop_type)*invRT) for qb in Qb[d][e] ]
            for d in range(i+4, j+1):
                Qs[i][j] += Qb[d][j]
            # Q recursion
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                if d == 0: # to deal with issue of wrapping around in the last iteration
                    Q[i][j] += Qs[d][j]
                else:
                    Q[i][j] += [ qs + q for qs in Qs[d][j] for q in Q[i][d-1] ]
    return sm.logsumexp(Q[0][N-1]) #np.sum(np.exp(Q[0][N-1]))

def linear_derivatives(param, sequence, N, g): # g : parameter that differentiating wrt
    return z.linear_derivatives_over_val(param, sequence, N, g)

def circular(param, sequence, N):
    g_stack = param[3]
    g_base_pair = gm.generator(sequence, param[0], param[1], param[2])

    # initializing general partition matrix
    Q = [[[0] for i in range(N)] for i in range(N)] # stores a,b,c, etc in exp(a) + exp(b) + ...
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    Qs = [[[] for i in range(N)] for i in range(N)]

    g_loop_over_RT = g_loop * invRT
    g_stack_over_RT = g_stack * invRT
    
    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                Qb[i][j].append(-(g_base_pair[i,j] + g_loop) * invRT)
            for d in xrange(i+1,j-4): # iterate over all possible rightmost pairs
                for e in xrange(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]: # interior loop possible
                            if i+1 == d and e+1 == j: #stacked
                                interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack
                            else: #loop
                                interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                            Qb[i][j] += [ qb - g_interior(g_base_pair[i,j], g_stack, interior_loop_type)*invRT for qb in Qb[d][e] ]
            for d in range(i+4, j+1):
                Qs[i][j] += Qb[d][j]
            # Q recursion
            if i == 0 and j == N-1: # closing chain
                for d in xrange(0, N-4):
                    if d == 0:
                        Q[i][j] += [ qs - g_loop_over_RT for qs in Qs[d][j] ]
                    else:
                        if len(Qb[0][d-1]) and len(Qb[d][N-1]): # to account for stacked pair forming when chain is closed
                            Q[i][j] += [ qs_k + q_m - g_loop_over_RT for qs_k in Qs[d][N-1] for q_m in Q[0][d-1] ]
                            for qs_k in Qs[d][N-1]:
                                for qs_w in Qs[0][d-1]:
                                    g_partial = qs_k + qs_w
                                    if (g_partial - g_loop*invRT) in Q[i][j]:
                                        Q[i][j].remove(g_partial - g_loop_over_RT) # removing single instance; acting as subtracting a term
                                        Q[i][j].append(g_partial - g_stack_over_RT)
                        else: # to account for interior loop forming when chain is closed
                            Q[i][j] += [ q + qs for q in Q[i][d-1] for qs in Qs[d][j] ]
            else:
                for d in xrange(i,j-3): # iterating over all possible rightmost pairs
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i][j] += Qs[d][j]
                    else:
                        Q[i][j] += [ qs + q for q in Q[i][d-1] for qs in Qs[d][j] ]
    return sm.logsumexp(Q[0][N-1])

def circular_derivatives(param, sequence, N, g):
    return z.circular_derivatives_over_val(param, sequence, N, g)
