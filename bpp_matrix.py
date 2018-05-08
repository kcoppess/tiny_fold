import tinyfold as tf
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

sns.set(font_scale=0.5)

ms2_hairpin = 'ACAUGAGGAUCACCCAUGU'

sequences = np.loadtxt('experimental/R101_Sequence.txt', dtype = 'string')
trouble = np.loadtxt('trouble.txt', dtype = 'string')
KD = np.loadtxt('experimental/R101_KDnoligand.txt', dtype = 'float')
#energies = np.array([-2.23026397, -4.61845988, -0.4459526 , -7.31024044, -4.43925842,
#       -3.63193668,  4.39840442, -1.52225387, -5.26723347,  1.34990893,
#       -5.72045913,  1.37062561])
energies = np.array([-4.00896773, -7.54098962, -0.23402838,  3.14820036, -5.44986465,
       -6.16460973,  2.81354172, -3.84467562, -6.34558658,  0.13278807,
       -4.91226052,  2.78108339])

for ii in range(100):
    SEQUENCE = sequences[ii]
    print ii
    start = SEQUENCE.index(ms2_hairpin)
    end = start + len(ms2_hairpin)

    tf_bpp = np.array(tf.RNA(SEQUENCE, False, list(energies), True, False).get_bpp_full())

    vienna_bpp = np.loadtxt("vienna_bpp_data/{}_bpp_matrix.txt".format(SEQUENCE))

    tf_bpp_matrix = pd.DataFrame(tf_bpp + tf_bpp.transpose(), index=list(SEQUENCE), columns=list(SEQUENCE))
    vienna_bpp_matrix = pd.DataFrame(vienna_bpp + vienna_bpp.transpose(), index=list(SEQUENCE), columns=list(SEQUENCE))
    
    subtitle_fs = 12.

    fig, ax = plt.subplots(figsize = (15, 5), ncols = 2, nrows = 1)
    if SEQUENCE in trouble:
        plt.suptitle("$KD_{exp}$ = %.4f"%KD[ii], x = 0.06, y = 0.5, fontsize = 10, color = 'r')
    else:
        plt.suptitle("$KD_{exp}$ = %.4f"%KD[ii], x = 0.06, y = 0.5, fontsize = 10, color = 'g')
    ax[0].xaxis.tick_top()
    ax[1].xaxis.tick_top()
    ax[0].set_title("tinyFold", fontsize = subtitle_fs)
    ax[1].set_title("Vienna", fontsize = subtitle_fs)
    ax[0].axhline(start, color = 'c', alpha=0.2)
    ax[1].axhline(start, color = 'c', alpha=0.2)
    ax[0].axhline(end, color = 'c', alpha=0.2)
    ax[1].axhline(end, color = 'c', alpha=0.2)
    ax[0].axvline(start, color = 'c', alpha=0.2)
    ax[1].axvline(start, color = 'c', alpha=0.2)
    ax[0].axvline(end, color = 'c', alpha=0.2)
    ax[1].axvline(end, color = 'c', alpha=0.2)

    sns.heatmap(tf_bpp_matrix, xticklabels=1, yticklabels=1, cmap="Greys", ax=ax[0], vmin = 0., vmax = 1.)
    sns.heatmap(vienna_bpp_matrix, xticklabels=1, yticklabels=1, cmap="Greys", ax=ax[1], vmin = 0., vmax = 1.)
    if SEQUENCE in trouble:
        plt.savefig('/Users/kcoppess/Desktop/tf_vienna/trouble/{}.png'.format(SEQUENCE))
    else:
        plt.savefig('/Users/kcoppess/Desktop/tf_vienna/fine/{}.png'.format(SEQUENCE))

#plt.show()
