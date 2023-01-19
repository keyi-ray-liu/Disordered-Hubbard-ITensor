import numpy as np
import matplotlib.pyplot as plt
import os


fig, ax = plt.subplots()

for f in os.listdir():

    if f[0] == 'e':

        val = float(f[6:])
        raw = np.loadtxt(f)

        ax.scatter( [val] * int(raw.shape[0]), raw )

plt.show()