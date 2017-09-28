# last modified: September 18, 2017
# calculates base-pair probability for particular base-pair in input RNA sequence
import numpy as np
import parameters as p

'''all free energy parameters in kcal/mol'''

h = p.h # differentiation step size

R = p.R # kcal/K/mol universal gas constant
T = p.T # K temperature (standard state - room temp)


# partition contribution from hairpin closed with base pair BP
def Q_hairpin(g_BP, g_loop):
    return np.exp(-(g_BP + g_loop)/(R*T))

# partition contribution from interior loop closed with base pair BP
# avoids double counting free energy from base pair formation
def Q_interior(g_BP, g_loop, g_stack, loop_type):
    if loop_type == 's':  # if stacked base pair
        return np.exp(-(g_BP + g_stack)/(R*T))
    elif loop_type == 'l':  # if loop
        return np.exp(-(g_BP + g_loop)/(R*T))
    else:
        print str(loop_type)+": Invalid type of interior loop"
        return 0.0

# required that base1 < base2
# generates matrix of base-pair probabilities
def mccaskill_linear(g_base_pair, g_loop, g_stack, N):
    # initializing general partition matrix
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0
    
    # initializing bound partition matrix
    Q_bound = np.zeros((N,N))
    
    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                Q_bound[i,j] = Q_hairpin(g_base_pair[i,j], g_loop)
            else: # no hairpin possible
                Q_bound[i,j] = 0.0
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            interior_loop_type = 's'
                        else: # if loop
                            interior_loop_type = 'l'
                        Q_bound[i,j] += Q_bound[d,e] * Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)
                    else: # no interior loop possible (one or both base pairs can't form)
                        Q_bound[i,j] += 0.0
            # Q recursion
            Q[i,j] = 1.0
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                for e in range(d+4,j+1):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i,j] += Q_bound[d,e]
                    else:
                        Q[i,j] += Q[i,d-1]*Q_bound[d,e]
    
    # base pair probability
    bp_prob = np.zeros((N,N))
    # need to start with outside pairs and work way in
    for i in range(N): # index for first base
        for j in range(N-1, -1, -1): # index for second base
            if i == 0 and j == N-1:
                bp_prob[i,j] = Q_bound[0, N-1] / Q[0, N-1]
            elif i == 0:
                bp_prob[i,j] = Q_bound[0, j] * Q[j+1, N-1] / Q[0, N-1]
            elif j == N-1:
                bp_prob[i,j] = Q[0, i-1] * Q_bound[i, N-1] / Q[0, N-1]
            else:
                bp_prob[i,j] = Q[0, i-1] * Q_bound[i,j] * Q[j+1, N-1] / Q[0, N-1] # if base-pair is not enclosed
                for k in range(i): #indexing outside bases
                    for l in range(j+1, N):
                        if k == i-1 and l == j+1:
                            interior_loop_type = 's'
                        else:
                            interior_loop_type = 'l'
                        if Q_bound[k,l]:
                            bp_prob[i,j] += bp_prob[k,l] * Q_bound[i,j] * Q_interior(g_base_pair[k,l], g_loop, g_stack, interior_loop_type) / Q_bound[k,l]
    return bp_prob

def flag_linear(base1, base2, g_base_pair, g_loop, g_stack, N):
    # initializing general partition matrix
    Q = [[[(1., False)] for i in range(N)] for i in range(N)] # stores exponential terms in partition function and whether the term is for a structure with base-pair in question
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    
    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                if i == base1 and j == base2:
                    dependent = True
                else:
                    dependent = False
                Qb[i][j].append((Q_hairpin(g_base_pair[i,j], g_loop), dependent))
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            interior_loop_type = 's'
                        else: # if loop
                            interior_loop_type = 'l'
                        for k in range(len(Qb[d][e])):
                            if Qb[d][e][k][1] or (i == base1 and j == base2) :
                                dependent = True
                            else:
                                dependent = False
                            Qb[i][j].append((Qb[d][e][k][0] * Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type), dependent))
            # Q recursion
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                for e in range(d+4,j+1):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i][j] += Qb[d][e]
                    else:
                        for k in range(len(Q[i][d-1])):
                            for m in range(len(Qb[d][e])):
                                if Q[i][d-1][k][1] or Qb[d][e][m][1]:
                                    dependent = True
                                else:
                                    dependent = False
                                Q[i][j].append((Q[i][d-1][k][0] * Qb[d][e][m][0], dependent))
    #print Q[0][N-1]
    part = 0.
    struc = 0.
    for f in Q[0][N-1]:
        part += f[0]
        if f[1]:
            struc += f[0]
    return struc/part

def flag_linear_derivatives(base1, base2, g_base_pair, g_loop, g_stack, N, g):
    # initializing general partition matrix
    Q = [[[(1., False)] for i in range(N)] for i in range(N)] # stores exponential terms in partition function and whether the term is for a structure with base-pair in question
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    
    # stores the derivatives
    dQ = [[[(0., False)] for i in range(N)] for i in range(N)] # stores derivative of corresponding term ^
    dQb = [[[] for i in range(N)] for i in range(N)]
    

    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j-i > 3 and g_base_pair[i,j]: # if possible hairpin: at least 4 positions apart and able to form a base pair
                if i == base1 and j == base2:
                    dependent = True
                else:
                    dependent = False
                Qb[i][j].append((Q_hairpin(g_base_pair[i,j], g_loop), dependent))
                if g == g_base_pair[i,j]:
                    dQb[i][j].append(((Q_hairpin(g_base_pair[i,j] + h, g_loop) - Q_hairpin(g_base_pair[i,j], g_loop))/h, dependent))
                else:
                    dQb[i][j].append((0., dependent))
            for d in range(i+1,j-4): # iterate over all possible rightmost pairs
                for e in range(d+4,j): # i < d < e < j and d,e must be at least 4 positions apart
                    interior_loop_type = ''
                    if g_base_pair[i,j] and g_base_pair[d,e]: # possible for both base pairs to form
                        if i+1 == d and e+1 == j: # if stacked
                            interior_loop_type = 's'
                        else: # if loop
                            interior_loop_type = 'l'
                        for k in range(len(Qb[d][e])):
                            if Qb[d][e][k][1] or (i == base1 and j == base2) :
                                dependent = True
                            else:
                                dependent = False
                            Qb[i][j].append((Qb[d][e][k][0] * Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type), dependent))
                            dQ_int = 0.
                            if g == g_base_pair[i,j]:
                                dQ_int = (Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                            elif g == g_stack:
                                dQ_int = (Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                            dQb[i][j].append(((dQ_int * Qb[d][e][k][0] + Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type) * dQb[d][e][k][0]), dependent))
            # Q recursion
            for d in range(i,j-3): # iterating over all possible rightmost pairs
                for e in range(d+4,j+1):
                    if d == 0: # to deal with issue of wrapping around in the last iteration
                        Q[i][j] += Qb[d][e]
                        dQ[i][j] += dQb[d][e]
                    else:
                        for k in range(len(Q[i][d-1])):
                            for m in range(len(Qb[d][e])):
                                if Q[i][d-1][k][1] or Qb[d][e][m][1]:
                                    dependent = True
                                else:
                                    dependent = False
                                Q[i][j].append((Q[i][d-1][k][0] * Qb[d][e][m][0], dependent))
                                dQ[i][j].append(((dQ[i][d-1][k][0] * Qb[d][e][m][0] + Q[i][d-1][k][0] * dQb[d][e][m][0]), dependent))
    part = 0.
    struc = 0.
    d_part = 0.
    d_struc = 0.
    for f in dQ[0][N-1]:
        d_part += f[0]
        if f[1]:
            d_struc += f[0]
    for f in Q[0][N-1]:
        part += f[0]
        if f[1]:
            struc += f[0]
    return (d_struc/part) - (struc * d_part / (part**2))


def mccaskill_circular(g_base_pair, g_loop, g_stack, N):
    Q = np.zeros((N,N))
    for m in range(N):
        Q[m,m-1] = 1.0

    Q_bound = np.zeros((N,N))

    for l in range(1, N+1): #length of subsequence
        for i in range(0, N-l+1): # start of subsequence
            j = i + l - 1
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                Q_bound[i,j] = Q_hairpin(g_base_pair[i,j], g_loop)
            else:  # hairpin not possible
                Q_bound[i,j] = 0.0
            for d in range(i+1,j-4):
                for e in range(d+4,j):
                    interior_loop_type = ''
                    if j - i > 3 and (i + N) - j > 3 and (d + N) - e > 3: #checking that the bases in each base pair are at least 4 positions apart
                        if g_base_pair[i,j] and g_base_pair[d,e]: # interior loop possible
                            if i+1 == d and e+1 == j: #stacked
                                interior_loop_type = 's' #g_interior = g_base_pair[i,j] + g_stack
                            else: #loop
                                interior_loop_type = 'l' #g_interior = g_base_pair[i,j] + g_loop
                            Q_bound[i,j] += Q_bound[d,e] * Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)
                        else: # interior loop not possible
                            Q_bound[i,j] += 0.0
                    else:
                        Q_bound[i,j] += 0.0
            Q[i,j] = 1.0
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    for e in range(d+4,N):
                        if d == 0:
                            Q[i,j] += Q_bound[0,e]*np.exp(-g_loop/(R*T))
                        else:
                            if e == N-1 and Q_bound[0,d-1] and Q_bound[d,N-1]: # to account for stacked pair forming when chain is closed
                                Q[i,j] += (Q[0,d-1] + Q_bound[0,d-1]*(np.exp(-(g_stack-g_loop)/(R*T)) - 1))*Q_bound[d,N-1]*np.exp(-g_loop/(R*T))
                            else: # to account for interior loop forming when chain is closed
                                Q[i,j] += Q[i,d-1]*Q_bound[d,e]*np.exp(-g_loop/(R*T))
            else:
                for d in range(i,j-3):
                    for e in range(d+4,j+1):
                        if d == 0:
                            Q[i,j] += Q_bound[d,e]
                        else:
                            Q[i,j] += Q[i,d-1]*Q_bound[d,e]
    
    # base pair probability
    bp_prob = np.zeros((N,N))
    # need to start with outside pairs and work way in
    for i in range(N): # index for first base
        for j in range(N-1, -1, -1): # index for second base
            prefactor = Q_bound[i,j] / Q[0, N-1]
            # term for if base pair is not enclosed
            bp_prob[i,j] = prefactor * np.exp(-g_loop/(R*T))
            for k in range(i): # left outside pairs
                for l in range(k+1, i):
                    if (k == j+1-N) and (l == i-1):
                        interior_loop_type = 's'
                        bp = prefactor * Q_bound[k,l] * np.exp(-g_stack/(R*T))
                    else:
                        interior_loop_type = 'l'
                        bp = prefactor * Q_bound[k,l] * np.exp(-g_loop/(R*T))
                    bp_prob[i,j] += bp
            for k in range(j+1, N): # right outside pairs
                for l in range(k+1, N):
                    if k == j+1 and l == i-1+N:
                        interior_loop_type = 's'
                        bp = prefactor * Q_bound[k,l] * np.exp(-g_stack/(R*T))
                    else:
                        interior_loop_type = 'l'
                        bp = prefactor * Q_bound[k,l] * np.exp(-g_loop/(R*T))
                    bp_prob[i,j] += bp 
            for k in range(i): #indexing outside bases
                for l in range(j+1, N):
                    if k == i-1 and l == j+1:
                        interior_loop_type = 's'
                    else:
                        interior_loop_type = 'l'
                    if Q_bound[k,l]:
                        bp_prob[i,j] += bp_prob[k,l] * Q_bound[i,j] * Q_interior(g_base_pair[k,l], g_loop, g_stack, interior_loop_type) / Q_bound[k,l]
    #print Q[0][N-1]
    return bp_prob


def flag_circular(base1, base2, g_base_pair, g_loop, g_stack, N):
    # initializing general partition matrix
    Q = [[[(1., False)] for i in range(N)] for i in range(N)] # stores a,b,c, etc in exp(a) + exp(b) + ...
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    
    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                if i == base1 and j == base2:
                    dependent = True
                else:
                    dependent = False
                Qb[i][j].append((Q_hairpin(g_base_pair[i,j], g_loop), dependent))
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
                                if Qb[d][e][k][1] or (i == base1 and j == base2):
                                    dependent = True
                                else:
                                    dependent = False
                                Qb[i][j].append(((Qb[d][e][k][0] * Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)), dependent))
            # Q recursion
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    for e in range(d+4,N):
                        if d == 0:
                            for k in range(len(Qb[d][e])):
                                if Qb[d][e][k][1]:
                                    dependent = True
                                else:
                                    dependent = False
                                Q[i][j].append(((Qb[d][e][k][0] * np.exp(-g_loop/(R*T))), dependent))
                        else:
                            if e == N-1 and len(Qb[0][d-1]) and len(Qb[d][N-1]): # to account for stacked pair forming when chain is closed
                                for k in range(len(Qb[d][N-1])):
                                    for m in range(len(Q[0][d-1])):
                                        if Qb[d][N-1][k][1] or Q[0][d-1][m][1]:
                                            dependent = True
                                        else:
                                            dependent = False
                                        Q[i][j].append(((Qb[d][N-1][k][0] * Q[0][d-1][m][0] * np.exp(-g_loop/(R*T))), dependent))
                                    for w in range(len(Qb[0][d-1])):
                                        if Qb[d][N-1][k][1] or Qb[0][d-1][w][1]:
                                            dependent = True
                                        else:
                                            dependent = False
                                        Q_partial = Qb[d][N-1][k][0] * Qb[0][d-1][w][0]
                                        Q_full = Q_partial * (np.exp(-g_stack/(R*T)) - np.exp(-g_loop/(R*T)))
                                        Q[i][j].append((Q_full, dependent))
                            else: # to account for interior loop forming when chain is closed
                                for k in range(len(Q[i][d-1])):
                                    for m in range(len(Qb[d][e])):
                                        if Q[i][d-1][k][1] or Qb[d][e][m][1]:
                                            dependent = True
                                        else:
                                            dependent = False
                                        Q[i][j].append((Q[i][d-1][k][0] * Qb[d][e][m][0] * np.exp(-g_loop/(R*T)), dependent))
            else:
                for d in range(i,j-3): # iterating over all possible rightmost pairs
                    for e in range(d+4,j+1):
                        if d == 0: # to deal with issue of wrapping around in the last iteration
                            Q[i][j] += Qb[d][e]
                        else:
                            for k in range(len(Q[i][d-1])):
                                for m in range(len(Qb[d][e])):
                                    if Q[i][d-1][k][1] or Qb[d][e][m][1]:
                                        dependent = True
                                    else:
                                        dependent = False
                                    Q[i][j].append((Q[i][d-1][k][0] * Qb[d][e][m][0], dependent))
    part = 0.
    struct = 0.
    for f in Q[0][N-1]:
        part += f[0]
        if f[1]:
            struct += f[0]
    #print struct/part
    return struct/part


def flag_circular_derivatives(base1, base2, g_base_pair, g_loop, g_stack, N, g):
    # initializing general partition matrix
    Q = [[[(1., False)] for i in range(N)] for i in range(N)] # stores a,b,c, etc in exp(a) + exp(b) + ...
    
    # initializing bound partition matrix
    Qb = [[[] for i in range(N)] for i in range(N)]
    
    # stores the derivatives
    dQ = [[[(0., False)] for i in range(N)] for i in range(N)] # stores derivative of corresponding term ^
    dQb = [[[] for i in range(N)] for i in range(N)]

    # calculation of partition function
    for l in range(1,N+1): # iterating over all subsequence lengths
        for i in range(0,N-l+1): # iterating over all starting positions for subsequences
            j = i + l - 1 # ending position for subsequence
            # Qb recursion
            if j - i > 3 and (i + N) - j > 3 and g_base_pair[i,j]: # checking that base pair can form and bases are at least 4 positions apart on both sides
                if i == base1 and j == base2:
                    dependent = True
                else:
                    dependent = False
                Qb[i][j].append((Q_hairpin(g_base_pair[i,j], g_loop), dependent))
                if g == g_base_pair[i,j]:
                    dQb[i][j].append(((Q_hairpin(g_base_pair[i,j] + h, g_loop) - Q_hairpin(g_base_pair[i,j], g_loop))/h, dependent))
                else:
                    dQb[i][j].append((0., dependent))
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
                                if Qb[d][e][k][1] or (i == base1 and j == base2):
                                    dependent = True
                                else:
                                    dependent = False
                                Qb[i][j].append(((Qb[d][e][k][0] * Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type)), dependent))
                                dQ_int = 0.
                                if g == g_base_pair[i,j]:
                                    dQ_int = (Q_interior(g_base_pair[i,j] + h, g_loop, g_stack, interior_loop_type) - Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                                elif g == g_stack:
                                    dQ_int = (Q_interior(g_base_pair[i,j], g_loop, g_stack + h, interior_loop_type) - Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type))/h
                                dQb[i][j].append(((dQ_int * Qb[d][e][k][0] + Q_interior(g_base_pair[i,j], g_loop, g_stack, interior_loop_type) * dQb[d][e][k][0]), dependent))
            # Q recursion
            if i == 0 and j == N-1: # closing chain
                for d in range(0, N-4):
                    for e in range(d+4,N):
                        if d == 0:
                            for k in range(len(Qb[d][e])):
                                if Qb[d][e][k][1]:
                                    dependent = True
                                else:
                                    dependent = False
                                Q[i][j].append(((Qb[d][e][k][0] * np.exp(-g_loop/(R*T))), dependent))
                                dQ[i][j].append((dQb[d][e][k][0] * np.exp(-g_loop/(R*T)), dependent))
                        else:
                            if e == N-1 and len(Qb[0][d-1]) and len(Qb[d][N-1]): # to account for stacked pair forming when chain is closed
                                for k in range(len(Qb[d][N-1])):
                                    for m in range(len(Q[0][d-1])):
                                        if Qb[d][N-1][k][1] or Q[0][d-1][m][1]:
                                            dependent = True
                                        else:
                                            dependent = False
                                        Q[i][j].append(((Qb[d][N-1][k][0] * Q[0][d-1][m][0] * np.exp(-g_loop/(R*T))), dependent))
                                        dQ[i][j].append(((dQb[d][N-1][k][0] * Q[0][d-1][m][0] + Qb[d][N-1][k][0] * dQ[0][d-1][m][0]) * np.exp(-g_loop/(R*T)), dependent))
                                    for w in range(len(Qb[0][d-1])):
                                        if Qb[d][N-1][k][1] or Qb[0][d-1][w][1]:
                                            dependent = True
                                        else:
                                            dependent = False
                                        Q_partial = Qb[d][N-1][k][0] * Qb[0][d-1][w][0]
                                        Q_full = Q_partial * (np.exp(-g_stack/(R*T)) - np.exp(-g_loop/(R*T)))
                                        Q[i][j].append((Q_full, dependent))
                                        #derivative:
                                        dQ_partial = dQb[d][N-1][k][0] * Qb[0][d-1][w][0] + Qb[d][N-1][k][0] * dQb[0][d-1][w][0]
                                        dQ_full = dQ_partial * (np.exp(-g_stack/(R*T)) - np.exp(-g_loop/(R*T)))
                                        if g == g_stack:
                                            dQ_full += Q_partial * np.exp(-g_stack/(R*T))/(-R*T)
                                        dQ[i][j].append((dQ_full, dependent))
                            else: # to account for interior loop forming when chain is closed
                                for k in range(len(Q[i][d-1])):
                                    for m in range(len(Qb[d][e])):
                                        if Q[i][d-1][k][1] or Qb[d][e][m][1]:
                                            dependent = True
                                        else:
                                            dependent = False
                                        Q[i][j].append((Q[i][d-1][k][0] * Qb[d][e][m][0] * np.exp(-g_loop/(R*T)), dependent))
                                        dQ[i][j].append(((dQ[i][d-1][k][0] * Qb[d][e][m][0] + Q[i][d-1][k][0] * dQb[d][e][m][0]) * np.exp(-g_loop/(R*T)), dependent))
            else:
                for d in range(i,j-3): # iterating over all possible rightmost pairs
                    for e in range(d+4,j+1):
                        if d == 0: # to deal with issue of wrapping around in the last iteration
                            Q[i][j] += Qb[d][e]
                            dQ[i][j] += dQb[d][e]
                        else:
                            for k in range(len(Q[i][d-1])):
                                for m in range(len(Qb[d][e])):
                                    if Q[i][d-1][k][1] or Qb[d][e][m][1]:
                                        dependent = True
                                    else:
                                        dependent = False
                                    Q[i][j].append((Q[i][d-1][k][0] * Qb[d][e][m][0], dependent))
                                    dQ[i][j].append(((dQ[i][d-1][k][0] * Qb[d][e][m][0] + Q[i][d-1][k][0] * dQb[d][e][m][0]), dependent))
    part = 0.
    struc = 0.
    d_part = 0.
    d_struc = 0.
    for f in dQ[0][N-1]:
        d_part += f[0]
        if f[1]:
            d_struc += f[0]
    for f in Q[0][N-1]:
        part += f[0]
        if f[1]:
            struc += f[0]
    return (d_struc/part) - (struc * d_part / (part**2))

