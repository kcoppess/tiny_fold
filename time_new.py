import time
import tinyfold as tf
import matplotlib.pyplot as plt
import numpy as np
import random

L = 50

sequence_length = np.arange(5,L+1)
n4_before = np.loadtxt("n4_partition_time_nonindex_gbp.txt", delimiter=',')
n4 = np.loadtxt("n4_partition_time.txt", delimiter=',')
n3 = np.loadtxt("n3_partition_time.txt", delimiter=',')

#plt.plot(sequence_length, n4_before, 'k', linewidth=2, label='N4 before')
plt.plot(sequence_length, n4, 'b',linewidth=2, label='N4')
plt.plot(sequence_length, n3, 'g',linewidth=2, label='N2')
plt.legend()
plt.xlabel("Sequence Length")
plt.ylabel("Run time (ms)")
plt.show()
