# last modified August 7, 2017
# version 1.0
# toy model to calculate partition function for a given sequence
# based off of NUPACK pseudo-code (N^4) given in Dirks & Pierce 2003
# TODO IGNORING MULTILOOPS
import numpy as np
import partition as z
import log_partition as log_z

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
        g_base_pair[m,n] = z.free_energy_pair(sequence[m],sequence[n],5.69,6.0,4.09)


if is_circular:
    print "Original: "+str(np.log(z.circular_sequence_partition(g_base_pair, g_loop, g_stack, N)))
    print "New vers: "+str(log_z.circular_sequence_logpartition(g_base_pair, g_loop, g_stack, N))
else:
    #print linear_sequence_partition(g_base_pair, N)
    print "Original: "+str(np.log(z.linear_sequence_partition(g_base_pair, g_loop, g_stack, N)))
    print "New vers: "+str(log_z.linear_sequence_logpartition(g_base_pair, g_loop, g_stack, N))
