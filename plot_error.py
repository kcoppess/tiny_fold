import matplotlib.pyplot as plt
import numpy as np

num_training_examples = 200

error = np.loadtxt("error.txt", delimiter=',')
size = np.arange(5, num_training_examples, 5)

plt.plot(size, error)
plt.show()
