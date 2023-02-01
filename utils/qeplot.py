import numpy as np
import matplotlib.pyplot as plt
import os
import sys


fig, ax = plt.subplots(figsize=(5, 10))
dirs = sys.argv[1]

files = os.listdir(dirs)
print(files)

for f in files:

    if f[0] == 'e':

        val = float(f[6:])
        raw = np.loadtxt(dirs + '/' + f)

        raw = raw - raw[0]

        ax.scatter( [val] * int(raw.shape[0]), raw, s=3 )

ax.set_xlabel('QE energy')
ax.set_ylabel('Energy level from GS')
ax.set_title('Plotting energy level as function of QE energy, ' + dirs)
fig.savefig( dirs + '.png')

plt.show()
