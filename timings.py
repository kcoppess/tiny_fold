import time
import partition as z
import base_pair_probabilities as bpp
import g_matrix as gm
import random
import matplotlib.pyplot as plt
import parameters as p

L = 50

sequence_length = range(5,L+1)
func_time = []
deri_time = []
grad_time = []
new_grad = []

num = 20

for l in range(5, L+1):
    start = time.clock()
    for i in range(20):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        g_base_pair = gm.generator(sequence, 5.69, 6., 4.09, l)
        part = bpp.mccaskill_linear(g_base_pair, 1., -7.09, l)
    end = time.clock()
    func_time.append((end-start)/10.)
    
    start = time.clock()
    for i in range(20):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        g_base_pair = gm.generator(sequence, 5.69, 6., 4.09, l)
        part = bpp.mccaskill_linear_derivatives(g_base_pair, 1., -7.09, l, 5.69)
    end = time.clock()
    deri_time.append((end-start)/10.)
    
    start = time.clock()
    for i in range(20):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        g_base_pair = gm.generator(sequence, 5.69, 6., 4.09, l)
        part = bpp.mccaskill_linear_derivatives(g_base_pair, 1., -7.09, l, 5.69)
        part = bpp.mccaskill_linear_derivatives(g_base_pair, 1., -7.09, l, 6.)
        part = bpp.mccaskill_linear_derivatives(g_base_pair, 1., -7.09, l, 4.09)
        part = bpp.mccaskill_linear_derivatives(g_base_pair, 1., -7.09, l, -7.09)
    end = time.clock()
    grad_time.append((end-start)/10.)

    #start = time.clock()
    #for i in range(20):
    #    sequence = ''.join(random.choice('AUGC') for _ in range(l))
    #    g_base_pair = gm.generator(sequence, 5.69, 6., 4.09, l)
    #    grad = grad_z.linear_derivatives(g_base_pair, p.energies,  1., l)
    #end = time.clock()
    #new_grad.append((end-start)/10.)

plt.plot(sequence_length, func_time, linewidth=2, label='Partition')
plt.plot(sequence_length, deri_time, linewidth=2, label='Single Derivative')
plt.plot(sequence_length, grad_time, linewidth=2, label='Full Gradient (Old)')
#plt.plot(sequence_length, new_grad, linewidth=2, label='Full Gradient (New)')
plt.xlabel('Sequence length')
plt.ylabel('Average Time (s)')
plt.legend(loc=2)
plt.show()
