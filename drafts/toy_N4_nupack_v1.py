# last modified August 7, 2017
# version 1.0
# toy model to calculate partition function for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# TODO IGNORING MULTILOOPS
import numpy as np

def free_energy_pair(b1,b2):
    g_AU = 5.69 # kcal/mol
    g_UG = 6.0  # kcal/mol
    g_GC = 4.09 # kcal/mol

    if (b1 == 'A' and b2 == 'U') or (b1 == 'U' and b2 == 'A'):
        return g_AU
    elif (b1 == 'U' and b2 == 'G') or (b1 == 'G' and b2 == 'U'):
        return g_UG
    elif (b1 == 'G' and b2 == 'C') or (b1 == 'C' and b2 == 'G'):
        return g_GC
    else:
        return 0.0

def linear_sequence_partition(g_base_pair, N):
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    
    # LINEAR SEQUENCE calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                g_hairpin = g_base_pair[i,j] + g_loop
                Qb[i,j] = np.exp(-g_hairpin/(R*T))
            else: # no hairpin possible
                g_hairpin = 0.0
                Qb[i,j] = 0.0
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    g_interior = 0.0
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            g_interior = g_base_pair[i,j] + g_stack # avoids double counting free energy from base pair formation
                        else: # if loop
                            g_interior = g_base_pair[i,j] + g_loop
                        Qb[i,j] += Qb[d,e] * np.exp(-g_interior/(R*T))
                    else: # no interior loop possible (one or both base pairs can't form)
                        Qb[i,j] += 0.0
            # Q recursion
            Q[i,j] = 1.0
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                for e in range(d+4, j+1):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i,j] += Qb[d,e]
                    else:
                        Q[i,j] += Q[i,d-1]*Qb[d,e]
    return Q[0,N-1]


# using finite difference method (basically analytical result)
def linear_sequence_derivatives(g_base_pair, N, g): # g : parameter that differentiating wrt
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    dQ = np.zeros((N,N)) # stores derivatives

    # initializing bound partition matrix
    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N)) # stores bound derivatives
    
    # LINEAR SEQUENCE calculation of partition function and its derivative
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                g_hairpin = g_base_pair[i,j] + g_loop
                Qb[i,j] = np.exp(-g_hairpin/(R*T))
                if g == g_base_pair[i,j] or g == g_loop:  # differentiating wrt current base-pair or loop parameter
                    d_g_hairpin = 1.0         # (g_hairpin + h - g_hairpin)/h
                    dQb[i,j] = -np.exp(-g_hairpin/(R*T)) * d_g_hairpin/(R*T)
                # else: dQb[i,j] = 0
            else: # no hairpin possible
                g_hairpin = 0.0
                Qb[i,j] = 0.0
                # dQb[i,j] = 0
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    g_interior = 0.0
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            g_interior = g_base_pair[i,j] + g_stack # avoids double counting free energy from base pair formation
                            if g == g_base_pair[i,j] or g == g_stack: # differentiating wrt base pair or stack parameter
                                d_g_interior = 1.0 # (g_interior + h - g_interior)/h
                                dQb[i,j] += -np.exp(-g_interior/(R*T)) * d_g_interior * Qb[d,e]/(R*T)
                        else: # if loop
                            g_interior = g_base_pair[i,j] + g_loop
                            if g == g_base_pair[i,j] or g == g_loop:
                                d_g_interior = 1.0
                                dQb[i,j] += -np.exp(-g_interior/(R*T)) * d_g_interior * Qb[d,e]/(R*T)
                        Qb[i,j] += Qb[d,e] * np.exp(-g_interior/(R*T))
                        dQb[i,j] += np.exp(-g_interior/(R*T)) * dQb[d,e]
                    else: # no interior loop possible (one or both base pairs can't form)
                        Qb[i,j] += 0.0
            # Q recursion
            Q[i,j] = 1.0
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                for e in range(d+4, j+1):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i,j] += Qb[d,e]
                        dQ[i,j] += dQb[d,e]
                    else:
                        Q[i,j] += Q[i,d-1]*Qb[d,e]
                        dQ[i,j] += dQ[i,d-1]*Q[d,e] + Q[i,d-1]*dQb[d,e]
    return Q[0,N-1], dQ[0,N-1]

def circular_sequence_partition(g_base_pair, N):
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0

    Qb = np.zeros((N,N))

    for l in range(1, N+1): #length of subsequence
        for i in range(0, N-l+1): # start of subsequence
            j = i + l - 1
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that bases are at least 4 positions apart on both sides
                g_hairpin = g_base_pair[i,j] + g_loop
                Qb[i,j] = np.exp(-g_hairpin/(R*T))
            else:
                g_hairpin = 0.0
                Qb[i,j] = 0.0
            for d in range(i+1,j-4):
                for e in range(d+4,j):
                    g_interior = 0.0
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]:
                            if i+1 == d and e+1 == j:
                                g_interior = g_base_pair[i,j] + g_stack
                            else:
                                g_interior = g_base_pair[i,j] + g_loop
                            Qb[i,j] += Qb[d,e]*np.exp(-g_interior/(R*T))
                        else:
                            Qb[i,j] += 0.0
                    else:
                        Qb[i,j] += 0.0
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    for e in range(d+4,N):
                        if d == 0:
                            Q[i,j] += Qb[0,e]*np.exp(-g_loop/(R*T))
                        else:
                            if e == N-1 and Qb[0,d-1] and Qb[d,N-1]: # to account for stacked pair forming when chain is closed
                                Q[i,j] += (Q[0,d-1] + Qb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) - 1))*Qb[d,N-1]*np.exp(-g_loop/(R*T))
                            else:
                                Q[i,j] += Q[i,d-1]*Qb[d,e]*np.exp(-g_loop/(R*T))
            else:
                for d in range(i,j-3):
                    for e in range(d+4,j+1):
                        if d == 0:
                            Q[i,j] += Qb[d,e]
                        else:
                            Q[i,j] += Q[i,d-1]*Qb[d,e]
    return Q[0,N-1]

sequence = raw_input('Sequence (if circular, end with \'-\'): ')
N = len(sequence)

is_circular = False
if sequence[-1] == '-':
    is_circular = True
    N = N - 1 # ignoring the '-' part of sequence


R = 0.0019872 # kcal/K/mol universal gas constant
T = 298.15 # K temperature (standard state - room temp)

# free energy parameters for loop closure and base pair stacking
g_loop = 1.0 # kcal/mol
g_stack = -7.09 # kcal/mol


# initializing base-pair free energy matrix
# IGNORING STERIC EFFECTS
g_base_pair = np.zeros((N,N))
for m in range(N):
    for n in range(N):
        g_base_pair[m,n] = free_energy_pair(sequence[m],sequence[n])


if is_circular:
    print circular_sequence_partition(g_base_pair, N)
else:
    #print linear_sequence_partition(g_base_pair, N)
    print linear_sequence_derivatives(g_base_pair, N,5.69)
