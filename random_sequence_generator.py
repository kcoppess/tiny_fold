# modified: August 9, 2017
# generates random RNA sequences and partition functions
import random
import partition as z
import numpy as np

f = open('sequences.txt','w')

training_data = []
for i in range(100):
    sequence = ''.join(random.choice('AUGC') for _ in range(random.randint(5,21)))
    N = len(sequence)
    datum = []
    datum.append(sequence)
    
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
        datum.append(z.circular_sequence_partition(g_base_pair, 1.0, -7.09, N)) # g_base_pair, g_loop, g_stack, N
    else:
        datum.append(z.linear_sequence_partition(g_base_pair, 1.0, -7.09, N)) # g_base_pair, g_loop, g_stack, N
    f.write(str(datum)+'\n')
    training_data.append(datum)

f.close()
