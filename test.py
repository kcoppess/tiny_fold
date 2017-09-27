import parameters as p
import numpy as np
import base_pair_probabilities as bpp
import partition as z
import matplotlib.pyplot as plt
import scipy.optimize as so
import g_matrix as gm

def func(param, sequence, b1, b2):
    N = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        N = N - 1 # ignoring the '-' part of sequence


    # free energy parameters for loop closure and base pair stacking
    g_loop = 1.0 # kcal/mol
    g_stack = param[3] # kcal/mol

    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    
    if is_circular:
        bpp_matrix = bpp.mccaskill_circular(g_base_pair, g_loop, g_stack, N)
    else:
        bpp_matrix = bpp.mccaskill_linear(g_base_pair, g_loop, g_stack, N)
    #print bpp_matrix[b1,b2]
    return bpp_matrix[b1, b2]

#sequence = 'AAACCCUUUCCC-'
#x = func([5.69, 6., 4.09, -7.09], sequence, 0, 3)

bpp1_au = []
bpp1_gu = []
bpp1_gc = []
bpp1_stack = []
bpp_flag_au = []
bpp_flag_gu = []
bpp_flag_gc = []
bpp_flag_stack = []

bpp_mccaskill = []
bpp_flag = []

num = 0

for sequence in p.training_sequences[5:]:
    #if num > 0:
    #    break

    M = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = True
        M = M - 1

    g_base_pair = gm.generator(sequence, 5.69, 6., 4.09, M)
    
    num += 1
    #print '-----------------------------------------------------------------------------'
    print sequence
    for i in range(M):
        for j in range(i+1, M):
            #print i, j
            bpp_mccaskill.append(func([5.69, 6., 4.09, -7.09], sequence, i, j))
            #grad = so.approx_fprime([5.69, 6., 4.09, -7.09], func, 1e-9, sequence, i, j)
            #bpp1_au.append(grad[0])
            #bpp1_gu.append(grad[1])
            #bpp1_gc.append(grad[2])
            #bpp1_stack.append(grad[3])
            if is_circular:
                f = bpp.flag_circular(i, j, g_base_pair, 1., -7.09, M)
            else:
                f = bpp.flag_linear(i, j, g_base_pair, 1., -7.09, M)
            bpp_flag.append(f)
            #bpp_flag.append(bpp.flag_circular(i, j, g_base_pair, 1., -7.09, M))
            #print f
            #print '.......'
            #bpp_flag_au.append(bpp.flag_linear_derivatives(i, j, g_base_pair, 1., -7.09, M, 5.69))
            #bpp_flag_gu.append(bpp.flag_linear_derivatives(i, j, g_base_pair, 1., -7.09, M, 6.))
            #bpp_flag_gc.append(bpp.flag_linear_derivatives(i, j, g_base_pair, 1., -7.09, M, 4.09))
            #bpp_flag_stack.append(bpp.flag_linear_derivatives(i, j, g_base_pair, 1., -7.09, M, -7.09))
            #y = bpp.flag_linear_derivatives(1, 7, g_base_pair, g_loop, g_stack, N, 4.09)
#x = np.linspace(min([min(bpp1_au), min(bpp1_gu),min(bpp1_gc),min(bpp1_stack)])-.5, max([max(bpp1_au), max(bpp1_gu),max(bpp1_gc),max(bpp1_stack)]) +.5, 100)
print np.linalg.norm(np.array(bpp_mccaskill) - np.array(bpp_flag))/np.sqrt(len(bpp_mccaskill))
x = np.linspace(min(bpp_mccaskill)-.5, max(bpp_flag)+.5, 100)

plt.scatter(bpp_mccaskill, bpp_flag)
#plt.scatter(bpp1_au, bpp_flag_au, color='blue', label='$g_{AU}$')
#plt.scatter(bpp1_gu, bpp_flag_gu, color='red', label='$g_{GU}$')
#plt.scatter(bpp1_gc, bpp_flag_gc, color='green', label='$g_{GC}$')
#plt.scatter(bpp1_stack, bpp_flag_stack, color='magenta', label='$g_{stack}$')
#plt.legend()
plt.plot(x, x, color='black')
plt.xlim([min(x), max(x)])
plt.ylim([min(x), max(x)])
plt.xlabel('McCaskill')
plt.ylabel('Boolean Flags')
plt.title('Base Pair Probabilities Two Ways (Linear/Circular)')
plt.show()
