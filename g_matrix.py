# generates free-energy matrix for an input sequence
# matrix used in several files
import numpy as np

# b1, b2: bases
# g_bp: free energy of forming that base pair (in kcal/mol)
def free_energy_pair(b1, b2, g_AU, g_UG, g_GC):
    if (b1 == 'A' and b2 == 'U') or (b1 == 'U' and b2 == 'A'):
        return g_AU
    elif (b1 == 'U' and b2 == 'G') or (b1 == 'G' and b2 == 'U'):
        return g_UG
    elif (b1 == 'G' and b2 == 'C') or (b1 == 'C' and b2 == 'G'):
        return g_GC
    else:
        return 0.0

def generator(sequence, g_au, g_ug, g_gc, N):
    g_base_pair = np.zeros((N,N))
    for n in range(N):
        for m in range(N):
            g_base_pair[n,m] = free_energy_pair(sequence[n], sequence[m], g_au, g_ug, g_gc)
    return g_base_pair
