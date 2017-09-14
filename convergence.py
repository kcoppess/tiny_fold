import matplotlib.pyplot as plt
import numpy as np
import log_partition as log_z
import stochastic_gradient_descent as sgd

K = sgd.K

adam_param = sgd.iteration_param

adam_Z = []

sequence = 'ACGCGGUGGGAUUCAC'
N = len(sequence)

energy_param = sgd.actual_param #[5.69, 6., 4.09, 1., -7.09]
gbp = np.zeros((N,N))
for m in range(N):
    for n in range(N):
        gbp[m,n] = log_z.free_energy_pair(sequence[m], sequence[n], energy_param[0], energy_param[1], energy_param[2])

actual = log_z.linear(gbp, energy_param[3], energy_param[4], N)
ACTUAL = []

for param in adam_param:
    g_base_pair = np.zeros((N,N))
    for m in range(N):
        for n in range(N):
            g_base_pair[m,n] = log_z.free_energy_pair(sequence[m], sequence[n], param[0], param[1], param[2])

    adam_Z.append((actual - log_z.linear(g_base_pair, param[3], param[4], N))**2 + 0.0001*np.dot(param, param))
    ACTUAL.append(sgd.a*np.dot(energy_param, energy_param))

plt.plot(K, adam_Z,linewidth = 2.0, label='Adam')
plt.plot(K, ACTUAL)
#plt.ylim(ymin = -500)
plt.title('Adaptive Learning Rate')
plt.xlabel('Iterations')
plt.ylabel('Cost Function')
plt.legend()
plt.show()
