import time
import partition as z
import base_pair_probabilities as bpp1
import random
import matplotlib.pyplot as plt
import parameters as p
import numpy as np

L = 50

sequence_length = range(5,L+1)
seq_len_3 = []
bpp = []
bpp_grad = []
part = []
grad = []

num = 20
'''
for l in range(5, L+1):
    seq_len_3.append(l**3/200.)
    start = time.clock()
    for i in range(num):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        parts = bpp1.mccaskill_linear(p.energies, sequence, l)#g_base_pair, 1., -7.09, l)
    end = time.clock()
    bpp.append(1e3*(end-start)/num)
    
    start = time.clock()
    for i in range(num):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        parts = bpp1.mccaskill_linear_gradient(p.energies, sequence, l) #g_base_pair, 1., -7.09, l, 5.69)
    end = time.clock()
    bpp_grad.append(1e3*(end-start)/num)
    
    start = time.clock()
    for i in range(num):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        parts = z.linear(p.energies, sequence, l) #g_base_pair, 1., -7.09, l, 5.69)
    end = time.clock()
    part.append(1e3*(end-start)/num)

    start = time.clock()
    for i in range(num):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        grads = z.linear_gradient(p.energies, sequence, l) #g_base_pair, p.energies,  1., l)
    end = time.clock()
    grad.append(1e3*(end-start)/num)
'''
c_bpp = np.loadtxt("C/bpp.txt", delimiter=' ', usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
c_bpp_grad = np.loadtxt("C/bpp_grad.txt", delimiter=' ', usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
c_part = np.loadtxt("C/part.txt", delimiter=' ', usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
c_grad = np.loadtxt("C/grad.txt", delimiter=' ', usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

cbpp = 1e-3*np.mean(c_bpp+c_part, axis=1)
cbpp_grad = 1e-3*np.mean(c_bpp_grad+c_part+c_grad+c_bpp, axis=1)
cpart = 1e-3*np.mean(c_part, axis=1)
cgrad = 1e-3*np.mean(c_grad+c_part, axis=1)

#plt.plot(sequence_length, seq_len_3, 'k')
plt.plot(sequence_length, cpart, 'b',linewidth=2, label='Partition')
plt.plot(sequence_length, cgrad, 'c',linewidth=2, label='Partition Gradient')
plt.plot(sequence_length, cbpp, 'r',linewidth=2, label='Basepair Probability (BPP)')
plt.plot(sequence_length, cbpp_grad, 'm',linewidth=2, label='BPP Gradient')
#plt.plot(sequence_length, part,'b--', linewidth=2)
#plt.plot(sequence_length, grad,'c--', linewidth=2)
#plt.plot(sequence_length, bpp,'r--', linewidth=2)
#plt.plot(sequence_length, bpp_grad,'m--', linewidth=2)
plt.xlabel('Sequence length')
plt.ylabel('Average Time (ms)')
plt.legend(loc=2)
plt.show()

