# last modified: August 16, 2017
# toy model to calculate NOTE log(partition function) for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# NOTE IGNORING MULTILOOPS
import numpy as np
import scipy.misc as sm
import partition as z

'''all free energy parameters in kcal/mol'''

h = 0.00000001 # differentiation step size

R = 0.0019872 # kcal/K/mol universal gas constant
T = 298.15 # K temperature (standard state - room temp)

# b1, b2: bases
# g_bp: free energy of forming that base pair in kcal/mol
def free_energy_pair(b1, b2, g_AU, g_UG, g_GC):
    if (b1 == 'A' and b2 == 'U') or (b1 == 'U' and b2 == 'A'):
        return g_AU
    elif (b1 == 'U' and b2 == 'G') or (b1 == 'G' and b2 == 'U'):
        return g_UG
    elif (b1 == 'G' and b2 == 'C') or (b1 == 'C' and b2 == 'G'):
        return g_GC
    else:
        return 0.0

# partition contribution from interior loop closed with base pair BP
# avoids double counting free energy from base pair formation
def g_interior(g_BP, g_loop, g_stack, loop_type):
    if loop_type == 's':  # if stacked base pair
        return g_BP + g_stack
    elif loop_type == 'l':  # if loop
        return g_BP + g_loop
    else:
        print str(loop_type)+": Invalid type of interior loop"
        return 0.0

def linear_sequence_partition(g_base_pair, g_loop, g_stack, N):
    # initializing general partition matrix
    Q = [[[0] for i in range(N)] for i in range(N)] # stores a,b,c, etc in exp(a) + exp(b) + ...
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    
    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                Qb[i][j].append(-(g_base_pair[i,j] + g_loop)/(R*T))
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            interior_loop_type = 's'
                        else: # if loop
                            interior_loop_type = 'l'
                        for k in range(len(Qb[d][e])):
                            Qb[i][j].append(Qb[d][e][k] - (g_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)/(R*T)))
            # Q recursion
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                for e in range(d+4,j+1):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i][j] += Qb[d][e]
                    else:
                        for k in range(len(Q[i][d-1])):
                            for m in range(len(Qb[d][e])):
                                Q[i][j].append(Q[i][d-1][k] + Qb[d][e][m])
    return sm.logsumexp(Q[0][N-1]) #np.sum(np.exp(Q[0][N-1]))


# using finite difference method
def linear_sequence_derivatives(g_base_pair, g_loop, g_stack, N, g): # g : parameter that differentiating wrt
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
                Qb[i,j] = z.Q_hairpin(g_base_pair[i,j], g_loop)
                if g == g_base_pair[i,j]:  # differentiating wrt current base pair parameter
                    dQb[i,j] = (z.Q_hairpin(g_base_pair[i,j] + h, g_loop) - z.Q_hairpin(g_base_pair[i,j], g_loop))/h
                elif g == g_loop: # differentiating wrt loop parameter
                    dQb[i,j] = (z.Q_hairpin(g_base_pair[i,j], g_loop + h) - z.Q_hairpin(g_base_pair[i,j], g_loop))/h
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
                        Qb[i,j] += Qb[d,e] * z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)
                        dQ_int = 0.0 # derivative of Q_interior
                        if g == g_base_pair[i,j]:
                            dQ_int = (z.Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                        elif g == g_loop and interior_loop_type == 'l':
                            dQ_int = (z.Q_interior(g_base_pair[i,j], g_loop + h, g_stack, interior_loop_type) - z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                        elif g == g_stack and interior_loop_type == 's':
                            dQ_int = (z.Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                        dQb[i,j] += dQ_int * Qb[d,e] + z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type) * dQb[d,e]
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
    return dQ[0,N-1]/z.linear_sequence_partition(g_base_pair, g_loop, g_stack, N)

def circular_sequence_partition(g_base_pair, g_loop, g_stack, N):
    # initializing general partition matrix
    Q = [[[0] for i in range(N)] for i in range(N)] # stores a,b,c, etc in exp(a) + exp(b) + ...
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    
    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                Qb[i][j].append(-(g_base_pair[i,j] + g_loop)/(R*T))
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]: # interior loop possible
                            if i+1 == d and e+1 == j: #stacked
                                interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack
                            else: #loop
                                interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                            for k in range(len(Qb[d][e])):
                                Qb[i][j].append(Qb[d][e][k] - (g_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)/(R*T)))
            # Q recursion
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    for e in range(d+4,N):
                        if d == 0:
                            for k in range(len(Qb[d][e])):
                                Q[i][j].append(Qb[d][e][k] - g_loop/(R*T))
                        else:
                            if e == N-1 and len(Qb[0][d-1]) and len(Qb[d][N-1]): # to account for stacked pair forming when chain is closed
                                for k in range(len(Qb[d][N-1])):
                                    for m in range(len(Q[0][d-1])):
                                        Q[i][j].append(Qb[d][N-1][k] + Q[0][d-1][m] - g_loop/(R*T))
                                    for w in range(len(Qb[0][d-1])):
                                        g_partial = Qb[d][N-1][k] + Qb[0][d-1][w]
                                        if (g_partial - g_loop/(R*T)) in Q[i][j]:
                                            Q[i][j].remove(g_partial - g_loop/(R*T)) # removing single instance; acting as subtracting a term
                                            Q[i][j].append(g_partial - g_stack/(R*T))
                            else: # to account for interior loop forming when chain is closed
                                for k in range(len(Q[i][d-1])):
                                    for m in range(len(Qb[d][e])):
                                        Q[i][j].append(Q[i][d-1][k] + Qb[d][e][m] - g_loop/(R*T))
            else:
                for d in range(i,j-3): # iterating over all possible rightmost pairs
                    for e in range(d+4,j+1):
                        if d == 0: # to deal with issue of wrapping around in the last iteration
                            Q[i][j] += Qb[d][e]
                        else:
                            for k in range(len(Q[i][d-1])):
                                for m in range(len(Qb[d][e])):
                                    Q[i][j].append(Q[i][d-1][k] + Qb[d][e][m])
    return sm.logsumexp(Q[0][N-1])


def circular_sequence_derivatives(g_base_pair, g_loop, g_stack, N, g):
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0

    dQ = np.zeros((N,N))

    Qb = np.zeros((N,N))
    dQb = np.zeros((N,N))

    for l in range(1, N+1): #length of subsequence
        for i in range(0, N-l+1): # start of subsequence
            j = i + l - 1
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                Qb[i,j] = z.Q_hairpin(g_base_pair[i,j], g_loop)
                if g == g_base_pair[i,j]:
                    dQb[i,j] = (z.Q_hairpin(g_base_pair[i,j] + h, g_loop) - z.Q_hairpin(g_base_pair[i,j], g_loop))/h
                elif g == g_loop:
                    dQb[i,j] = (z.Q_hairpin(g_base_pair[i,j], g_loop + h) - z.Q_hairpin(g_base_pair[i,j], g_loop))/h
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
                            Qb[i,j] += Qb[d,e] * z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)
                            dQ_int = 0.0
                            if g == g_base_pair[i,j]:
                                dQ_int = (z.Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                            elif g == g_loop and interior_loop_type == 'l':
                                dQ_int = (z.Q_interior(g_base_pair[i,j], g_loop + h, g_stack, interior_loop_type) - z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                            elif g == g_stack and interior_loop_type == 's':
                                dQ_int = (z.Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                            dQb[i,j] += dQ_int * Qb[d,e] + z.Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type) * dQb[d,e]
                        else: # interior loop not possible
                            Qb[i,j] += 0.0
                    else:
                        Qb[i,j] += 0.0
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    for e in range(d+4,N):
                        if d == 0:
                            Q[i,j] += Qb[0,e]*np.exp(-g_loop/(R*T))
                            if g == g_loop:
                                dQ[i,j] += dQb[0,e] * np.exp(-g_loop/(R*T)) - Qb[0,e] * np.exp(-g_loop/(R*T))/(R*T)
                            else:
                                dQ[i,j] += dQb[0,e] * np.exp(-g_loop/(R*T))
                        else:
                            if e == N-1 and Qb[0,d-1] and Qb[d,N-1]: # to account for stacked pair forming when chain is closed
                                Q[i,j] += (Q[0,d-1] + Qb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) - 1))*Qb[d,N-1]*np.exp(-g_loop/(R*T))
                                if g == g_loop:
                                    dQ[i,j] += (dQ[0,d-1] + dQb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) -1) + Qb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T))/(R*T)))*Qb[d,N-1]*np.exp(-g_loop/(R*T))
                                    dQ[i,j] += (Q[0,d-1] + Qb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) - 1))*(dQb[d,N-1]*np.exp(-g_loop/(R*T)) - Qb[d,N-1]*np.exp(-g_loop/(R*T))/(R*T))
                                elif g == g_stack:
                                    dQ[i,j] += (dQ[0,d-1] + dQb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) -1) + Qb[0,d-1]*(-np.exp(-(g_stack-g_loop)/(R*T))/(R*T)))*Qb[d,N-1]*np.exp(-g_loop/(R*T))
                                    dQ[i,j] += (Q[0,d-1] + Qb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) - 1))*dQb[d,N-1]*np.exp(-g_loop/(R*T))
                                else:
                                    dQ[i,j] += (dQ[0,d-1] + dQb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) -1))*Qb[d,N-1]*np.exp(-g_loop/(R*T))
                                    dQ[i,j] += (Q[0,d-1] + Qb[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) - 1))*dQb[d,N-1]*np.exp(-g_loop/(R*T))
                            else: # to account for interior loop forming when chain is closed
                                Q[i,j] += Q[i,d-1]*Qb[d,e]*np.exp(-g_loop/(R*T))
                                if g == g_loop:
                                    dQ[i,j] += dQ[i,d-1]*Qb[d,e]*np.exp(-g_loop/(R*T)) + Q[i,d-1]*dQb[d,e]*np.exp(-g_loop/(R*T)) - Q[i,d-1]*Qb[d,e]*np.exp(-g_loop/(R*T))/(R*T)
                                else:
                                    dQ[i,j] += dQ[i,d-1]*Qb[d,e]*np.exp(-g_loop/(R*T)) + Q[i,d-1]*dQb[d,e]*np.exp(-g_loop/(R*T))
            else:
                for d in range(i,j-3):
                    for e in range(d+4,j+1):
                        if d == 0:
                            Q[i,j] += Qb[d,e]
                            dQ[i,j] += dQb[d,e]
                        else:
                            Q[i,j] += Q[i,d-1]*Qb[d,e]
                            dQ[i,j] += dQ[i,d-1] * Qb[d,e] + Q[i,d-1] * dQb[d,e]
    return dQ[0,N-1]/z.circular_sequence_partition(g_base_pair, g_loop, g_stack, N)