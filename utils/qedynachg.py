import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from scipy.fftpack import fft, ifft

def compplot( timestep=[0.1]):

    def animate(i):

        print("frame {}".format(i))
        prev = -1
        cnt = 0
        colors = ['blue', 'purple', 'magenta', 'green']
        colorsQElim = ['black', 'brown', 'grey']
        colorsQE = ['red', 'orange', 'yellow']

        for j in range(num):
            
            # internal counter in case multiple data on one plot
            if ax_pos[j] != prev:
                cnt =0 
                prev = ax_pos[j]
                cur_title = inputs[j]
                ax[0][ax_pos[j]].clear()
                ax[1][ax_pos[j]].clear()

            else:
                cnt += 1
                cur_title = cur_title + ' and ' +  inputs[j]

            #ax[j].set_ylim(lo[j], hi[j])

            # set QE axis
            ax[0][ax_pos[j]].set_ylim(0, 1)
            #ax[j].scatter( [1, L[j] + 4 ], data[j][i, [0, L[j] + 3]], c='red', label=r'$QE \   \ |0\rangle $')

            # plot QE in separated axis
            qe_label = [ 'Left QE', 'Right QE' ]
            ticks_loc = [0 , 1]
            ax[0][ax_pos[j]].set_title(cur_title + 'time = {:2f}'.format(timestep[j] * i) + 'avg {}'.format(avg[j]))
            ax[0][ax_pos[j]].scatter( [0, 1], data[j][i, [1, L[j] + 2]], c=colorsQE[cnt])

            ax[0][ax_pos[j]].scatter( [0], leftlo[j], c=colorsQElim[cnt], label='maxQE_{}'.format(inputs[j]))
            ax[0][ax_pos[j]].scatter( [1], righthi[j], c=colorsQElim[cnt], label='maxQE_{}'.format(inputs[j]))

            ax[0][ax_pos[j]].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            ax[0][ax_pos[j]].set_xticklabels(qe_label)

            x = np.arange(L[j] + avg[j] )
            y = data[j][i, 2:L[j]+2] - gs[j]

            # adding 0 paddings on both ends

            if avg[j]:
                padding = np.zeros(1)
                y = (np.concatenate( (padding, y), axis=0) + np.concatenate((y, padding) , axis=0) ) / 2

            ax[1][ax_pos[j]].scatter( x, y, c=colors[cnt], label='Chain_{}'.format(inputs[j])) 
            ax[1][ax_pos[j]].plot(x, y, c=colors[cnt])
            ax[1][ax_pos[j]].set_xlabel('sites')
            ax[1][ax_pos[j]].set_ylim(lo[j], hi[j])

            ax[1][ax_pos[j]].set_ylabel( r'$ \langle n_i \rangle $')
            ax[1][ax_pos[j]].xaxis.set_major_locator(MaxNLocator(integer=True))

        #ax[1][-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        fig.tight_layout()
        #fig.suptitle('time = {:2f}'.format(timestep * i), color=color)


    # num of input cases

    inputs = np.loadtxt('input_files', dtype=str)
    timestep = np.loadtxt('timestep')
    avg = np.loadtxt('avg', dtype=int)
    ax_pos = np.loadtxt('ax', dtype=int)

    num = inputs.shape[0]
    outdir = os.getcwd() + '/vids/' 

    if num != avg.shape[0] or num!= timestep.shape[0] or num!= ax_pos.shape[0]:
        raise ValueError('input dim does not match')
    
    # update: separated QE and actual
    fig, ax = plt.subplots(2, ax_pos[-1] + 1, figsize=(5 * num + 1,9))

    if ax_pos[-1] + 1 == 1:
        ax = [ [ax[0]],[ax[1]]]

    dirs = [''] * num
    data = [[]] * num
    gs = [[]] * num
    frames = [0] * num
    L = [0] * num
    hi = np.zeros(num)
    lo = np.zeros(num)
    
    leftlo = np.ones(num)
    righthi = np.zeros(num)
    
    for i in range(num):
        dirs[i] = os.getcwd() + '/' + inputs[i ] + '/'
        data[i] = np.loadtxt(dirs[i] + 'expN')

        int_tag = 'nonint' if 'nonint' in inputs[i] else ''

        frames[i], L[i] = data[i].shape
        L[i] -= 4

        gs[i] = np.loadtxt( os.getcwd() + '/gs{}{}'.format(L[i], int_tag))

        padding = np.zeros( (frames[i], 1))
        temp_y = data[i][:, 2:L[i]+ 2] -  np.tile(gs[i], (frames[i], 1))

        if avg[i]:
            temp_y = (np.concatenate( (padding, temp_y), axis=1) + np.concatenate((temp_y, padding) , axis=1) ) / 2

        hi[i] = np.amax(temp_y)
        lo[i] = np.amin(temp_y)

        leftlo[i] = np.amin( data[i][:, 1])
        righthi[i] = np.amax( data[i][:, L[i] + 2])

    #box = ax[-1].get_position()
    #ax[-1].set_position([box.x0, box.y0, box.width * 0.9, box.height])

    frame = min(frames)

    print(frame)
    anim = FuncAnimation(fig, animate, frames= frame, interval=100, repeat=False)
    #plt.show()

    #mpl.rcParams['animation.ffmpeg_path'] = os.getcwd() + '/ffmpeg'
    signature = [ str(x) for x in list(inputs) + ['avg'] + list(avg)]

    html = False
    if not html:
        writervideo = animation.FFMpegWriter(fps=15)
        anim.save( outdir + 'Dyna{}.mp4'.format( '_'.join(signature)), writer=writervideo)

    else:
        writervideo = animation.HTMLWriter()
        anim.save( outdir + 'Dyna{}.html'.format( '_'.join(signature)), writer=writervideo)

def compplot2d( timestep=0.1, timechange=1e10, comment=""):

    def animate(i):

        for j in range(num):
            
            for k in range(3):
                ax[j * 3 + k].clear()

            cax[j].clear()

            raw = data[j][i]

            chain = raw [ 2:-2]
            chain = np.reshape(chain, shapes[j])

            ax[j * 3].scatter( [0], [raw[0]], label='QE1gs')
            ax[j * 3].scatter( [1], [raw[1]], label='QE1ex')
            ax[j * 3].set_ylim(0, 1)

            ax[j * 3].set_xticks([0, 1], ['QE1gs', 'QE1ex'])

            ax[j * 3 + 2].scatter( [1], [raw[-1]], label='QE2gs')
            ax[j * 3 + 2].scatter( [0], [raw[-2]], label='QE2ex')

            ax[j * 3 + 2].set_xticks([0, 1], ['QE2ex', 'QE2gs'])
            ax[j * 3 + 2].set_ylim(0, 1)
            

            im = ax[j *3 + 1].imshow( chain)
            
            if i * timestep < timechange:
                color = 'black'

            else:
                color = 'red'

            ax[j * 3 + 1].set_title('time = {:2f}'.format(timestep * i), color=color)
            ax[j * 3 + 1].axis('off')
            fig.colorbar(im, cax=cax[j], orientation='horizontal', label=r'$\langle n \rangle$')
            
            

        fig.suptitle(comment)

    input = sys.argv[2:]

    outdir = os.getcwd() + '/vids/' 
    num = len(input) 
    fig, ax = plt.subplots( num , 3, figsize=(15, 5 * num + 1), width_ratios=[1, 10, 1] * num)


    dirs = [''] * num
    data = [[]] * num
    frames = [0] * num
    shapes = []
    cax = []
    
    for i in range(num):
        dirs[i] = os.getcwd() + '/' + input[i] + '/'
        data[i] = np.loadtxt(dirs[i] + 'expN')

        frames[i], shape = data[i].shape
        shape = (2, shape//2 - 2)


        divider = make_axes_locatable(ax[i * 3 + 1])
        cax.append( divider.append_axes('bottom', size='10%',  pad=0.05) )

        
        shapes.append( shape)

    print(cax)
    print(shapes)

    box = ax[-1].get_position()
    ax[-1].set_position([box.x0, box.y0, box.width * 0.9, box.height])

    frame = min(frames)

    print(frame)
    anim = FuncAnimation(fig, animate, frames= frame, interval=100, repeat=False)
    #plt.show()

    #mpl.rcParams['animation.ffmpeg_path'] = os.getcwd() + '/ffmpeg'

    writervideo = animation.FFMpegWriter(fps=15)
    anim.save( outdir + 'Dyna{}.mp4'.format( '_'.join(input)), writer=writervideo)

    #writervideo = animation.HTMLWriter()
    #anim.save( outdir + 'Dyna{}.html'.format( '_'.join(input[1:])), writer=writervideo)


def bonds():

    def func(x, a, b):
        return a * np.exp(b* x) + bond[0]
    
    file = sys.argv[2]
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


def ft():

    file_dir = os.getcwd() + '/' + sys.argv[2] + '/'
    data= np.loadtxt(file_dir + 'expN')

    tseries = data[:,3]
    X = fft(tseries)

    #print(X)
    N = len(X)
    allfreq = np.arange(N)
    posfreq = np.arange(1, N//2 + 1)

    fig, ax = plt.subplots(4, figsize=(10, 20))

    ax[0].scatter(np.arange(tseries.shape[0]), tseries)
    ax[0].set_title('original oscialltion vs. t')
    ax[0].set_ylabel( r'$ \langle \hat{n}_i \rangle $')
    ax[0].set_xlabel( 'time (1/hopping)')

    ax[1].stem(posfreq, np.abs(X[1: N //2 + 1]), 'b', \
            markerfmt=" ", basefmt="-b")
    ax[1].set_xlabel('Freq (1/t)')
    ax[1].set_ylabel('FFT Amplitude')
    ax[1].set_title('FFT')


    ax[2].plot( allfreq, ifft(X))
    ax[2].set_title('Inverse FFT')
    ax[2].set_ylabel( r'$ \langle \hat{n}_i \rangle $')
    ax[2].set_xlabel( 'time (1/hopping)')

    ax[3].scatter( 1/posfreq * N, np.abs(X[1: N//2 + 1]))
    ax[3].set_ylabel('FFT Amplitude')
    ax[3].set_xlabel( 'time interval of pulse (1/t)')

    fig.tight_layout()
    fig.savefig('fft.pdf')


def overlap():

    file = sys.argv[2]

    len = file[:2]
    overlap = np.loadtxt('overlap/' + file)

    fig, ax = plt.subplots()
    im = ax.imshow( np.log10(overlap.transpose()))
    fig.colorbar(im, label=r'$log_{10} \ |\langle eigen| \psi \rangle| $')

    ax.set_xlabel('time')
    ax.set_ylabel('eigen state- time independent')
    ax.set_title('Overlap between QE eigenstates and TE state from GS: {} system'.format(len))
    plt.show()

def overlap_gs():

    file = sys.argv[2]

    overlap = np.loadtxt( file + '/overlap')

    fig, ax = plt.subplots()
    ax.scatter( np.arange(overlap.shape[0]), overlap)

    ax.set_xlabel('Eigenstate')
    ax.set_ylabel('Overlap')
    ax.set_title('Overlap between QE eigenstates and TE_GS')
    
    fig.savefig( file + '/overlap.png')


if __name__ == '__main__':
    
    control = int(sys.argv[1])

    timechange = 1e10
    #comment = "Timechange at {}, Left: bond dim capped at 150. Right: bond dim capped at 300".format(timechange)

    if control == 1:
        compplot( )

    elif control == 2:
        compplot2d(timestep=0.1)

    elif control == 3:
        bonds()

    elif control == 4:
        overlap()

    elif control == 5:
        ft()

    elif control == 6:
        overlap_gs()