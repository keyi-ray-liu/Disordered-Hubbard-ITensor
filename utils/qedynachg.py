import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit

def compplot( timestep=0.1, timechange=1e10, comment=""):

    def animate(i):

        for j in range(num):
            ax[j].clear()
            ax[j].set_ylim(lo[j], hi[j])
            ax[j].scatter( [1, L[j] + 4 ], data[j][i, [0, L[j] + 3]], c='red', label=r'$QE \   \ |0\rangle $')
            ax[j].scatter( [2, L[j] + 3], data[j][i, [1, L[j] + 2]], c='green', label=r'$QE\   \ |1\rangle $')

            ax[j].scatter( np.arange(3, 3 + L[j] ), data[j][i, 2:L[j]+2], c='blue', label='Chain')

            ax[j].set_xlabel('sites')

            if i * timestep < timechange:
                color = 'black'

            else:
                color = 'red'

            ax[j].set_title('time = {:2f}'.format(timestep * i), color=color)
            ax[j].set_ylabel( r'$ \langle n_i \rangle $')
            ax[j].xaxis.set_major_locator(MaxNLocator(integer=True))

        ax[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        fig.suptitle(comment)

    input = sys.argv

    outdir = os.getcwd() + '/vids/' 
    num = len(input) -1 
    fig, ax = plt.subplots(1, num, figsize=(5 * num + 1,5))

    if num == 1:
        ax = [ax]

    dirs = [''] * num
    data = [[]] * num
    frames = [0] * num
    L = [0] * num
    hi = [0] * num
    lo = [0] * num
    
    for i in range(num):
        dirs[i] = os.getcwd() + '/' + input[i+1] + '/'
        data[i] = np.loadtxt(dirs[i] + 'expN')

        frames[i], L[i] = data[i].shape
        L[i] -= 4

        hi[i] = np.amax(data[i])
        lo[i] = np.amin(data[i])

    box = ax[-1].get_position()
    ax[-1].set_position([box.x0, box.y0, box.width * 0.9, box.height])

    frame = min(frames)

    print(frame)
    anim = FuncAnimation(fig, animate, frames= frame, interval=100, repeat=False)
    #plt.show()

    #mpl.rcParams['animation.ffmpeg_path'] = os.getcwd() + '/ffmpeg'

    writervideo = animation.FFMpegWriter(fps=15)
    anim.save( outdir + 'Dyna{}.mp4'.format( '_'.join(input[1:])), writer=writervideo)

    #writervideo = animation.HTMLWriter()
    #anim.save( outdir + 'Dyna{}.html'.format( '_'.join(input[1:])), writer=writervideo)


def bonds():

    def func(x, a, b):
        return a * np.exp(b* x) + bond[0]
    
    file = sys.argv[1]
    bond = np.loadtxt('bonds/' + file)

    time = np.linspace(0, (len(bond) - 1 )/10, len(bond))

    print(time)


    fig, ax = plt.subplots()
    ax.scatter(time, bond, label='Data')
    ax.set_title('Bond dimension vs. time,' + file)
    ax.set_xlabel('Time')
    ax.set_ylabel('Max. bond dim. of the MPS')

    popt, pcov = curve_fit(func, time, bond)

    ax.scatter(time, func(time, *popt), label='{:.2f} * exp ( {:.2f} * t) + C(t=0)'.format(*popt))

    ax.legend()

    plt.show()


def overlap():

    file = sys.argv[1]

    len = file[:2]
    overlap = np.loadtxt('overlap/' + file)

    fig, ax = plt.subplots()
    im = ax.imshow( np.log10(overlap.transpose()))
    fig.colorbar(im, label=r'$log_{10} \ |\langle eigen| \psi \rangle| $')

    ax.set_xlabel('time')
    ax.set_ylabel('eigen state- time independent')
    ax.set_title('Overlap between QE eigenstates and TE state from GS: {} system'.format(len))
    plt.show()


if __name__ == '__main__':

    timechange = 13.4
    #comment = "Timechange at {}, Left: bond dim capped at 150. Right: bond dim capped at 300".format(timechange)
    comment = "cutoff small"
    compplot(timestep=0.1, timechange=timechange, comment=comment)
    #bonds()
    #overlap()