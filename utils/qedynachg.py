import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib import axes
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from scipy.fftpack import fft, ifft, fftfreq

def compplot( dim= 1 ):

    def animate(i):

        print("frame {}".format(i))
        prev = -1
        cnt = 0
        colors = ['blue', 'red', 'green', 'orange']
        colorsQElim = ['black', 'brown', 'grey']
        colorsQE = ['red', 'orange', 'yellow']

        for j in range(num):
            
            # internal counter in case multiple data on one plot
            if ax_pos[j] != prev or dim >1:

                cnt =0 
                prev = ax_pos[j]
                cur_title = inputs[j]
                ax[0][ax_pos[j]].clear()
                ax[1][ax_pos[j]].clear()

                if dim > 1:
                    cax[j].clear()

            else:
                cnt += 1
                cur_title = cur_title + ' & ' +  inputs[j]

            #ax[j].set_ylim(lo[j], hi[j])
            frame_cnt = (i + 1) * skip_frame[j] - 1
            # set QE axis
            ax[0][ax_pos[j]].set_ylim(0, 1)
            #ax[j].scatter( [1, L[j] + 4 ], data[j][i, [0, L[j] + 3]], c='red', label=r'$QE \   \ |0\rangle $')

            # plot QE in separated axis
            qe_label = [ 'Left QE', 'Right QE' ]
            ticks_loc = [0 , 1]
            ax[0][ax_pos[j]].set_title(cur_title + 'time = {:.2f}'.format(timestep[j] * (frame_cnt + 1)) + 'avg {}'.format(avg[j]))
            ax[0][ax_pos[j]].scatter( [0, 1], data[j][frame_cnt, [1, L[j] + 2]], c=colorsQE[cnt])

            ax[0][ax_pos[j]].scatter( [0], leftlo[j], c=colorsQElim[cnt], label='min_left_QE_{}'.format(inputs[j]))
            ax[0][ax_pos[j]].scatter( [1], righthi[j], c=colorsQElim[cnt], label='max_right_QE_{}'.format(inputs[j]))

            ax[0][ax_pos[j]].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            ax[0][ax_pos[j]].set_xticklabels(qe_label)

            y = data[j][frame_cnt, 2:L[j]+2] - gs[j]
            if dim == 1:
                x = np.arange(L[j] + avg[j] )
                

                # adding 0 paddings on both ends

                if avg[j]:
                    padding = np.zeros(1)
                    y = (np.concatenate( (padding, y), axis=0) + np.concatenate((y, padding) , axis=0) ) / 2

                ax[1][ax_pos[j]].scatter( x, y, c=colors[cnt], label='Chain_{}'.format(inputs[j]), s=3) 
                ax[1][ax_pos[j]].plot(x, y, c=colors[cnt])
                ax[1][ax_pos[j]].set_xlabel('sites')
                ax[1][ax_pos[j]].set_ylim(lo[j], hi[j])

                ax[1][ax_pos[j]].set_ylabel( r'$ \langle n_i \rangle $')
                ax[1][ax_pos[j]].xaxis.set_major_locator(MaxNLocator(integer=True))

                ax[0][ax_pos[j]].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=1)
                ax[1][ax_pos[j]].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=1)

            else:

                side = int(np.sqrt( y.shape[0]))
                y = np.reshape(y, (side, side))
                im = ax[1][ax_pos[j]].imshow( y, vmin=lo[j], vmax=hi[j])
                
                #ax[1][ax_pos[j]].axis('off')
                fig.colorbar(im, cax=cax[j], orientation='horizontal', label=r'$\langle n \rangle$')


        #ax[1][-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        fig.tight_layout()
        #fig.suptitle('time = {:2f}'.format(timestep * i), color=color)


    # num of input cases

    inputs = np.loadtxt('input_files', dtype=str, ndmin=1)
    timestep = np.loadtxt('timestep', ndmin=1)
    avg = np.loadtxt('avg', dtype=int, ndmin=1)
    ax_pos = np.loadtxt('ax', dtype=int, ndmin=1)
    skip_frame = np.loadtxt( 'skip_frame', dtype=int, ndmin=1)
    gs_tag = np.loadtxt('gs_tag', dtype=str, ndmin=1)

    num = inputs.shape[0]

    outdir = os.getcwd() + '/vids/' 

    if num != avg.shape[0] or num!= timestep.shape[0] or num!= ax_pos.shape[0] or num!= skip_frame.shape[0] or num!= gs_tag.shape[0]:
        raise ValueError('input dim does not match')
    
    # update: separated QE and actual
    #fig, ax = plt.subplots(2, ax_pos[-1] + 1, figsize=(4 * num + 1, 12))
    fig, ax = plt.subplots(2, ax_pos[-1] + 1, figsize=(32, 18))

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

    if dim >1:
        cax = []
    
    for i in range(num):
        dirs[i] = os.getcwd() + '/' + inputs[i ] + '/'
        data[i] = np.loadtxt(dirs[i] + 'expN')

        frames[i], L[i] = data[i].shape
        L[i] -= 4

    common_frame = min( [ frames[i] // skip_frame[i] for i in range(num)])
    print(common_frame)

    for i in range(num):
        gs[i] = np.loadtxt( os.getcwd() + '/gs{}'.format(gs_tag[i]))

        padding = np.zeros( (frames[i], 1))
        temp_y = data[i][:, 2:L[i]+ 2] -  np.tile(gs[i], (frames[i], 1))

        if avg[i]:
            temp_y = (np.concatenate( (padding, temp_y), axis=1) + np.concatenate((temp_y, padding) , axis=1) ) / 2

        hi[i] = np.amax(temp_y)
        lo[i] = np.amin(temp_y)

        cutoff = min(common_frame * skip_frame[i], frames[i])
        leftlo[i] = np.amin( data[i][:cutoff, 1])
        righthi[i] = np.amax( data[i][:cutoff, L[i] + 2])

        if dim >1:
            divider = make_axes_locatable(ax[1][ax_pos[i]])
            cax.append( divider.append_axes('bottom', size='10%',  pad=0.05) )

    #box = ax[-1].get_position()
    #ax[-1].set_position([box.x0, box.y0, box.width * 0.9, box.height])


    anim = FuncAnimation(fig, animate, frames= common_frame, interval=100, repeat=False)
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


def expecteigen():

    file = sys.argv[2]
    raw = np.loadtxt( file + '/expeigen')
    
    row = 3
    col = 5

    fig, ax = plt.subplots(row, col, figsize = (5 * col, 5 * row))

    cnt = 0
    for ex in raw:

        r = cnt // col
        c = cnt % col

        cur_ax : axes.Axes = ax[r][c]

        chain = np.arange(1, len(ex) - 3)
        qe = [0, len(ex) - 3]
        cur_ax.scatter(chain, ex[2:-2], color = 'red', s= 5)
        cur_ax.scatter( qe, ex[ [1, -1] ], color ='orange')

        cur_ax.set_ylabel(r'$ \langle \hat{n} \rangle$')
        cur_ax.set_xlabel('site')

        cur_ax.set_title('Eigen state {}'.format(cnt))

        cnt += 1
    fig.savefig( 'plots/eigen{}.pdf'.format(file))


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

    def dephasing(t, a, c, d, off):
        return a * np.cos(  t * c) * np.exp( - t/d) + off
    
    
    file_dir = os.getcwd() + '/' + sys.argv[2] + '/'
    plot_dir = os.getcwd() + '/vids/'
    data= np.loadtxt(file_dir + 'expN')

    timestep = float(sys.argv[3])
    T = 1/timestep
    tseries = data[:,1]
    N = tseries.shape[0]
    
    times = np.arange(N) * timestep

    
    X = fft(tseries)
    paras, _ = curve_fit(dephasing, times, tseries, [1, 0.01, 1000, 0.5])
    #print(X)
    N = len(X)
    #posfreq = np.arange(1, N//2 + 1)

    posfreq = fftfreq(N, T)[1:N//2+1]
    

    fig, ax = plt.subplots(3, figsize=(10, 15))

    ax[0].scatter(times, tseries)
    ax[0].plot( times, dephasing(times, *paras), label='${:.2f} \cos ({:.5f} * t) * \exp (-t/{:.2f}) + {:.2f}$'.format(*paras), color='red')
    ax[0].set_title('original oscialltion vs. t')
    ax[0].set_ylabel( r'$ \langle \hat{n}_i \rangle $')
    ax[0].set_xlabel( 'time (1/hopping)')
    ax[0].legend()

    #ax[1].stem(posfreq, np.abs(X[1: N //2 + 1]), 'b', markerfmt=" ", basefmt="-b")

    ax[1].plot(posfreq, np.abs(X[1: N //2 + 1]))
    ax[1].set_xlim(0, 0.2)
    ax[1].set_xlabel('Freq (1/t)')
    ax[1].set_ylabel('FFT Amplitude')
    ax[1].set_title('FFT')


    ax[2].scatter( 1/posfreq * N, np.abs(X[1: N//2 + 1]))
    ax[2].set_ylabel('FFT Amplitude')
    ax[2].set_xlabel( 'time interval of pulse (1/t)')

    fig.tight_layout()
    fig.savefig(plot_dir + 'fft_{}.pdf'.format(sys.argv[2]))


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
        compplot( dim=1)

    elif control == 2:
        compplot( dim =2)

    elif control == 3:
        bonds()

    elif control == 4:
        overlap()

    elif control == 5:
        ft()

    elif control == 6:
        overlap_gs()

    elif control == 7:
        expecteigen()