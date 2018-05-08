import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

bpp = np.loadtxt("bpp_vienna.txt")

sequence = "GAAAGGGCAGGAUAUGCAGAGAGAGAGAGAGAGAGAAAAAAACAUGAGGAUCACCCAUGUAGAAGGGCCGAAAG"

n = len(sequence)

bpp_matrix = np.zeros((n,n))

for entry in bpp:
    bpp_matrix[int(entry[0])-1, int(entry[1])-1] = entry[2]

bpp_m = pd.DataFrame(bpp_matrix, index=list(sequence), columns=list(sequence))
ax = sns.heatmap(bpp_matrix, xticklabels=1, yticklabels=1, cmap="Greys")
ax.xaxis.tick_top()
plt.show()

