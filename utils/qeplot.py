import numpy as np
import matplotlib.pyplot as plt
import os
import sys


fig, ax = plt.subplots()
dirs = sys.argv[1]

files = os.listdir(dirs)
print(files)

for f in files:

    if f[0] == 'e':

        val = float(f[6:])
        raw = np.loadtxt(dirs + '/' + f)

        raw = raw - raw[0]

        ax.scatter( [val] * int(raw.shape[0]), raw, s=1 )

plt.show()
