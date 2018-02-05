import tinyfold as tf
import filtering as fil
import numpy as np

#energies = [ 6.33147038,  6.03501064,  5.21450578, -5.69452531] # w = 0.000001
#energies = [ 8.00662263,  6.07959333,  5.61716355, -5.45409029] # w = 1e-11
#energies = [  4.09736972, -23.39842545, -15.42124886,  16.69237574] # w = 1e-20, tol = 1e-15
#energies = [ -2.31331946, -13.93911931, -20.90039423,  13.23790764] # new tol = 1e-18
energies = [ -2.35309106, -14.09522637, -20.16892953,  13.63350101] # new tol = 1e-20

'''BPP TRAINING DATA'''
num_training_examples = 60
rna = []
actual_bpp = []
guess = []
#energy_param = p.energies
sequences = fil.Sequence[:num_training_examples] #p.training_sequences[:10]
for i in range(num_training_examples):
    rna.append(tf.RNA(sequences[i], False, list(energies), True, False))
    bp = fil.closing_bp_indices[i]
    guess.append(rna[i].get_bpp(bp[0], bp[1]))
    actual_bpp.append((1e-9) / fil.KDnoligand[i])
print 'finished gathering training data'

print np.linalg.norm((np.array(actual_bpp[:50]) - np.array(guess[:50])) / np.array(actual_bpp[:50]))
print np.linalg.norm((np.array(actual_bpp[50:]) - np.array(guess[50:])) / np.array(actual_bpp[50:]))
print np.array(np.array(actual_bpp) - np.array(guess))
