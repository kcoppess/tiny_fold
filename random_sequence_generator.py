# generates random RNA sequences
import random

N = 5

f = open('sequences.txt','w')

for i in range(10):
    sequence = ''.join(random.choice('AUGC') for _ in range(N))
    f.write(str(sequence)+'\n')

f.close()
