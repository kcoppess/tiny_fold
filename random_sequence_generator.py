# modified: August 27, 2017
# generates random RNA sequences and log(partition)
import random
import log_partition as log_z
import numpy as np

# num_train: number of training sequences to be generated
# energy_param: [g_AU, g_GU, g_GC, g_loop, g_stack] in kcal/mol
# minlen: minimum length of random sequences
# maxlen: maximum length of random sequences
def training_generator(num_train, energy_param, minlen, maxlen):
    f = open('sequences_train.txt','w')

    training_data = np.zeros(num_train) # storing values of partition function
    sequences = []              # corresponding RNA sequence
    for i in range(num_train):
        sequence = ''.join(random.choice('AUGC') for _ in range(random.randint(minlen, maxlen + 1)))
        if random.randint(0,1000) < 500: #generating circular sequences
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
                g_base_pair[m,n] = log_z.free_energy_pair(sequence[m], sequence[n], energy_param[0], energy_param[1], energy_param[2])
        
        if is_circular:
            training_data[i] = log_z.circular(g_base_pair, energy_param[3], energy_param[4], N)
        else:
            training_data[i] = log_z.linear(g_base_pair, energy_param[3], energy_param[4], N)
        f.write(str(sequence)+' '+str(training_data[i])+'\n')
    
    f.close()
    return training_data, sequences

# num_test: number of test sequences to be generated
# energy_param: [g_AU, g_GU, g_GC, g_loop, g_stack] in kcal/mol
# minlen: minimum length of random sequences
# maxlen: maximum length of random sequences
def test_generator(num_test, energy_param, minlen, maxlen):
    f = open('sequences_test.txt','w')

    test_data = np.zeros(num_test) # storing values of partition function
    sequences = []              # corresponding RNA sequence
    for i in range(num_test):
        sequence = ''.join(random.choice('AUGC') for _ in range(random.randint(minlen, maxlen + 1)))
        if random.randint(0,1000) < 500: #generating circular sequences
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
                g_base_pair[m,n] = log_z.free_energy_pair(sequence[m], sequence[n], energy_param[0], energy_param[1], energy_param[2])
        
        if is_circular:
            test_data[i] = log_z.circular(g_base_pair, energy_param[3], energy_param[4], N)
        else:
            test_data[i] = log_z.linear(g_base_pair, energy_param[3], energy_param[4], N)
        f.write(str(sequence)+' '+str(test_data[i])+'\n')
    
    f.close()
    return test_data, sequences
