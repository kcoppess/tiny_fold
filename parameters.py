# holds parameter values for model
import numpy as np

'''-------FREE ENERGY PARAMETERS-----'''
# free energy for forming loop
g_loop = 1.

# energy parameters used to generate the training/test data
# [g_AU, g_GU, g_GC, g_stack]
energies = np.array([5.69, 6., 4.09, -7.09])

'''-------COST FUNCTION PARAMETERS-----'''
# weight on prior term
w = 1.

# priors
priors = np.array([5., 6., 4., -7])

'''--------LEARNING PARAMETERS------'''
# Adam parameters
alpha = 0.0001
beta1 = 0.9
beta2 = 0.999

'''-----------DATA SETS------------'''
training_data = np.array([  2.03865532e-01,   1.74345206e-03,   1.48952522e-02,   4.01105700e+00,  4.53177624e-02,   8.15228162e-04,   5.75364405e-03,   1.33490248e+00, 2.27751039e-02,   1.79234946e-02,   2.00795684e-03,   3.97134027e-02, 2.55239052e-03,   2.55648910e-02,   6.46582613e-02,   6.01486963e-05,  1.76533659e+00,   6.52938822e-01,   3.25810173e+00,   1.30279590e+00,  1.71729096e+00,   2.50493160e-03,   2.83040490e-02,   2.30034625e-03,  1.79665217e-02,   8.41719002e-01,   4.13533559e-01,   4.11844464e-02,  7.29529870e+00,   6.08295156e-04,   4.04648593e-02,   8.68634121e-04,  1.05216623e-03,   6.52232625e-02,   6.38816606e-01,   6.30687827e-02,  7.49505711e-03,   1.07769642e-02,   1.88855038e-02,   9.77691159e-02,  2.18154125e+00,   4.78485599e-02,   1.08599479e-02,   5.76867661e-04,  5.69087713e-03,   2.61834734e-02,   4.30795410e+00,   3.33627984e+00,  5.57212365e+00,   7.46901161e-05])

training_sequences = ['CCAAAUGUCGACCGGCUUCAAU', 'CCAAACUAAAAUCUAAAGAGC-', 'AACACCGAGUCUAACCCCAUGA', 'GUUAAUGGUUGGGGGGCUUCCC-', 'CUCAAUUUCCAUAGCUUAUAGG', 'AUCUAUCGAAG-', 'CCUAGCAACACAUAUCAGC-', 'UAUAUACGUAUCCGGUAUA', 'UGGACAACACUCAGG', 'AGCUAACAUAAGGGGU', 'AAGGGAUUCG', 'AGCGGACACCUCAAAG', 'GUACACCACA', 'ACCAAGAUCCUUC', 'AGCGCUGAUCUGAGCU-', 'GAGGAUAAGGGGAU-', 'UCUGCGCUCUCCGCCUGUAUC', 'CCCCGCUAACGGCUCA-', 'CACUUACCGGAGAGUGGCUA', 'GGGGGGUGCUACGGGUACGG-', 'UUAACCUCGUGGUUGUUC-', 'GUCGGAACAUA', 'CUCAUUCAUUGAUUUU', 'CAAGGCAUUUAUUAA', 'UGUGAGAUUUAUUGCGUG-', 'CGGAGGGUUCGCUUGGCGUUAU-', 'GGUAACCUGCAUAAAAUAGG-', 'UGCUCAGAGGCCCCAC-', 'AGGCACCCCCUGCCGCCGC-', 'UUUACGACUCUGGCAAA-', 'UGAGUCAUUGAUAUG', 'GUGUGAAAAUG', 'CAGGUCAUCACGCA-', 'UUCCGAGCGUCAGGAGU-', 'GCAGCGCUCUCGGUUACC-', 'CAAACAUGCAUUUUGC-', 'CUAUGAGUUUAUUUAAAUCCG-', 'CGUCGCGCAAAUUGU-', 'UUAGAUCUGAGUACUUUGAUAA-', 'GACCCGGCUUUGAAAGGUGUAC-', 'AUAUGACUCCUAGUCACAUAAG-', 'UGUCGAUAGUCCGUC', 'CCACGACGUUGACUUCUAA', 'CCGUUGUCACC-', 'UGCAUCCAGUGAAGA-', 'CUUAAGCUGCUUCCGA', 'CAUUGGUCCUGGCGUCGCUGC-', 'UCGAUGACAUCGGUGGUCUCU', 'UUUCAUGUAAACAUGAG', 'CUACACUCACG-']

test_data = np.array([  1.11408639e-03,   1.66383002e-03,   5.01944001e+00,   5.19297796e-03,  6.28788222e+00])

test_sequences = ['GCACGCAAGGAA', 'UAUGGUGCGGU', 'CUGCCCCCGUCAGCGGGUC-', 'CCCUAAUAACUAGAAGA-', 'GGUAGUGCUCUACCGCAACGG']