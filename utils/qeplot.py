import numpy as np
import matplotlib.pyplot as plt
import os
import sys


fig, ax = plt.subplots()
dirs = sys.argv[1]

print(os.listdir(dirs))

for f in os.listdir():

    if f[0] == 'e':

        val = float(f[6:])
        raw = np.loadtxt(f)

        raw = raw - raw[0]

        ax.scatter( [val] * int(raw.shape[0]), raw )

plt.show()