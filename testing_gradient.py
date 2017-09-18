import numpy as np
import parameters as p
import matplotlib.pyplot as plt
import train_toy as train
import scipy.optimize as so

random_vectors = []
random_magnitudes = []
current = []

error = []
approx = []

g = np.linspace(-7.5,-2.0,40)

for i in g:
    randvec = np.zeros(4) #random.uniform(-8, 0, size=(4))
    randvec[0] = 5.69
    randvec[1] = 6.
    randvec[2] = 4.09
    randvec[3] = i
    random_vectors.append(randvec)
    random_magnitudes.append(randvec[3])

    #app = so.approx_fprime(randvec, train.cost, 1e-8, p.g_loop, p.training_data, p.training_sequences)
    #approx.append(np.linalg.norm(app))
    
    curr = train.cost_gradient(randvec, p.g_loop, p.training_data, p.training_sequences)
    current.append(curr[3])

    #err = so.check_grad(train.cost, train.cost_gradient, randvec, p.g_loop, p.training_data, p.training_sequences)
    #error.append(err)
    #print str(randvec) + '  ' + str(np.linalg.norm(randvec)) + '  ' + str(err)

#x = range(int(max(current)))

plt.plot(random_magnitudes, current)
#plt.plot(x, x)
plt.xlabel('Random Parameter')
plt.ylabel('Gradient')
plt.show()
