import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.widgets import Slider

def gen(site, par):

    if par == 0:

        return [[0] * site]

    if site == 0:

        if par > 0:
            return []

        return [[]]


    cur = []

    for comb in gen(site - 1, par - 1):
        cur.append( [1] + comb)

    for comb in gen(site -1, par):
        cur.append( [0] + comb)

    return cur

def draw(fig, ax, site, par, conf):


    width = 0.1
    offset = 0.1
    

    col = site
    row = par * 2 + 1
    
    for i in range(0, row, 2):
        ax.scatter( np.arange(col), np.ones(col) * i, facecolors='none', edgecolors='black',s=60)

    for j in range(1, row, 2):
        ax.scatter( np.arange(col), np.ones(col) * j, c='black', s=60)

    for arr in conf:
        
        cur = arr[0]

        #print(arr)
        for j in range(len(arr) -1 ):

            if arr[j] != arr[j + 1]:

                ax.arrow( j + offset, cur , 1 - 3 * offset , 1 - 2 * offset, head_width=width)
                cur += 1

            elif arr[j] == 0:
                ax.arrow( j + offset, cur, 1 - 3 * offset, 0, head_width=width)

            else:
                ax.arrow( j+ offset, cur, 1 - 3 * offset, 2 - 2 * offset, head_width=width)
                cur += 2
            
    ax.invert_yaxis()
    ax.axis('off')
    plt.draw()
    plt.show()
    

def main():
    ax1 = plt.axes([0.25, 0.1 ,0.5, 0.02])
    ax2 = plt.axes([0.25, 0.2 ,0.5, 0.02])
    sites = Slider(ax1, 'No. of site', 1, 12, valstep=1)
    pars = Slider(ax2, 'No. of fermions', 1, 12, valstep=1)

    def update(val):
        ax.clear()
        site = int(sites.val)
        par = int(pars.val)

        conf = gen(site, par)
        draw(fig, ax, site, par, conf)

    fig, ax = plt.subplots()
    sites.on_changed(update)
    pars.on_changed(update)
    plt.show()

main()






