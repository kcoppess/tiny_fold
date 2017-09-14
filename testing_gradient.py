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

for i in range(50):
    randvec = np.random.uniform(-10, 10, size=(4))
    random_vectors.append(randvec)
    random_magnitudes.append(np.linalg.norm(randvec))

    #app = so.approx_fprime(randvec, train.cost, 1e-8, p.g_loop, p.training_data, p.training_sequences)
    #approx.append(np.linalg.norm(app))
    
    #curr = train.cost_gradient(randvec, p.g_loop, p.training_data, p.training_sequences)
    #current.append(np.linalg.norm(curr))

    err = so.check_grad(train.cost, train.cost_gradient, randvec, p.g_loop, p.training_data, p.training_sequences)
    error.append(err)
    print str(randvec) + '  ' + str(np.linalg.norm(randvec)) + '  ' + str(err)

#x = range(int(max(current)))

plt.scatter(random_magnitudes, error)
#plt.plot(x, x)
plt.xlabel('Magnitude of Parameter Vector')
plt.ylabel('Error for log(Z)-gradient')
plt.show()
