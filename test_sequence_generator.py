# modified: August 14, 2017
# generates random RNA sequences and partition functions to test trained model on
import random
import partition as z
import numpy as np

f = open('sequences_test.txt','w')

M = 5

g_AU = 5.69 # kcal/mol
g_GU = 6.0  # kcal/mol
g_GC = 4.09 # kcal/mol
g_loop = 1.0 # kcal/mol
g_stack = -7.09 # kcal/mol

test_data = np.zeros(M) # storing values of partition function
sequences = []              # corresponding RNA sequence
for i in range(M):
    sequence = ''.join(random.choice('AUGC') for _ in range(random.randint(5,21)))
    if random.randint(0,1000) < 500: # circular sequences
        sequence = sequence + '-'
    N = len(sequence)
    sequences.append(sequence)
    
    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring the '-' part of sequence
        
    # initializing base-pair free energy matrix
    # IGNORING STERIC EFFECTS
    g_base_pair = np.zeros((N,N))                                                    
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = z.free_energy_pair(sequence[m], sequence[n], g_AU, g_GU, g_GC)
                                                                             
    if is_circular:
        test_data[i] = z.circular_sequence_partition(g_base_pair, g_loop, g_stack, N)
    else:
        test_data[i] = z.linear_sequence_partition(g_base_pair, g_loop, g_stack, N)
    f.write(str(sequence)+' '+str(test_data[i])+'\n')

f.close()
