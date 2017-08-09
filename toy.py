# last modified: August 8, 2017
# version 3.0
# toy model to calculate partition function for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# TODO IGNORING MULTILOOPS
import numpy as np
import partition as z
'''all free energy parameters in kcal/mol'''

f = open('sequences.txt','r')

p = open('training_data.txt','w')

sequence = f.readline()[:-1]
N = len(sequence)
while sequence:
    p.write(str(sequence))

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring the '-' part of sequence
    
    # initializing base-pair free energy matrix
    # IGNORING STERIC EFFECTS
    g_base_pair = np.zeros((N,N))
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = z.free_energy_pair(sequence[m],sequence[n],5.69,6.0,4.09) # bp1, bp2, g_AU, g_GU, g_GC
    
    if is_circular:
        p.write(' '+str(z.circular_sequence_partition(g_base_pair, 1.0, -7.09, N))) # g_base_pair, g_loop, g_stack, N
        for dg in [5.69, 6.0, 4.09, 1.0, -7.09]:
            p.write(' '+str(z.circular_sequence_derivatives(g_base_pair, 1.0, -7.09, N, dg)))
    else:
        p.write(' '+str(z.linear_sequence_partition(g_base_pair, 1.0, -7.09, N))) # g_base_pair, g_loop, g_stack, N
        for dg in [5.69, 6.0, 4.09, 1.0, -7.09]:
            p.write(' '+str(z.linear_sequence_derivatives(g_base_pair, 1.0, -7.09, N, dg))) # g_base_pair, g_loop, g_stack, N, g
    
    p.write('\n')
    sequence = f.readline()[:-1]
    N = len(sequence)

f.close()
p.close()
