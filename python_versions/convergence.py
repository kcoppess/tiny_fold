# produces graph of value of cost function after each iteration of SGD
import matplotlib.pyplot as plt
import numpy as np
#import log_partition as log_z
import stochastic_gradient_descent as sgd
import parameters as p
import base_pair_probabilities as bpp

K = sgd.K

adam_param = sgd.iteration_param

adam_Z = []

sequence = 'ACGCGGUGGGAUUCAC'
N = len(sequence)

energy_param = p.energies #[5.69, 6., 4.09, 1., -7.09]

actual = bpp.mccaskill_linear(energy_param, sequence, N)[6,10]
ACTUAL = []

for param in adam_param:
    prior = p.priors - param
    adam_Z.append((actual - bpp.mccaskill_linear(param, sequence, N)[6,10])**2 + p.w*np.dot(prior, prior))
    ACTUAL.append(p.w*np.dot(p.priors - energy_param, p.priors - energy_param))

plt.plot(K, adam_Z,linewidth = 2.0, label='Adam')
plt.plot(K, ACTUAL)
#plt.ylim(ymin = -500)
plt.title('Convergence')
plt.xlabel('Iterations')
plt.ylabel('Cost Function')
plt.legend()
plt.show()
