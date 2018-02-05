import time
import tinyfold as tf
import matplotlib.pyplot as plt
import numpy as np
import random

L = 50

sequence_length = np.arange(5,L+1)
seq_len_3 = []
bpp = []
bpp_grad = []
part = []
grad = []

num = 100

energy = [5.69, 6., 4.09, -7.09]

for l in range(5, L+1):
    seq_len_3.append(7.6e-5*l**2 + 0.01)
    start = time.clock()
    for i in range(num):
        sequence = ''.join(random.choice('AUGC') for _ in range(l))
        parts = tf.RNA(sequence, False, energy, False, False)
    end = time.clock()
    bpp.append(1e3*(end-start)/num)

np.savetxt("n3_partition_time.txt", bpp, delimiter=',')

plt.plot(sequence_length, seq_len_3, 'k')
plt.plot(sequence_length, bpp, 'b',linewidth=2)
plt.show()
