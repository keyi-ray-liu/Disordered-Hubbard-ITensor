import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from scipy.fft import fft, rfft, fftfreq, fftshift
import glob
from matplotlib import colormaps
from matplotlib.widgets import Button, Slider

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
        gs[i] = np.loadtxt( os.getcwd() + '/{}'.format(gs_tag[i]))

        if len(gs[i]) != L[i]:
            gs[i] = gs[i][2:-2]

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

        cur_ax : plt.Axes = ax[r][c]

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

    files = sys.argv[2:]

    fig, ax = plt.subplots()

    for file in files:
        overlap = np.loadtxt( file + '/overlap')
        overlap = np.abs(overlap.flatten()) ** 2

        ax.scatter( np.arange(overlap.shape[0]), overlap, label=file[-3:])
        ax.plot( np.arange(overlap.shape[0]), overlap, label=file[-3:])

        ax.set_xlabel('Eigenstate')
        ax.set_ylabel('$ |overlap|^2$')
        ax.set_title('Overlap between QE eigenstates and TE_GS')
    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
    fig.savefig( 'plots/overlap{}.png'.format('_'.join(files)))

def check_pl_level():

    file = sys.argv[2]

    fig, ax = plt.subplots()

    ex = np.loadtxt( file )
    ex = ex.flatten()

    ex = (ex - ex[0]) / (ex[1] - ex[0])

    ax.scatter( np.arange(ex.shape[0]), ex)
    ax.plot( np.arange(ex.shape[0]), ex)

    ax.set_xlabel('Eigenstate')
    ax.set_ylabel('Pl num')
    ax.set_title('Plasmon number vs. eigenstate')
    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.savefig( 'plots/check_pl{}.png'.format(file) )

def plot_corr():

    def animate(i):
        print("frame {}".format(i))

        # axes has shape (r, c) where each col represents a case
        for j, row in enumerate(axes):
            for k, ax in enumerate(row):
                
                #print(j, k)
                cax : plt.Axes = caxes[j][k]
                ax.clear()
                cax.clear()

                im = ax.imshow( Corrs[j, k, i], vmin = lo[j,k], vmax=  hi[j,k])
                
                #ax[1][ax_pos[j]].axis('off')
                ax.set_title( vals[j] + files[k])
                fig.colorbar(im, cax=cax,  orientation='horizontal' )



    files = sys.argv[3:]
    outdir = os.getcwd() + '/vids/' 
    
    op = sys.argv[2]

    if op == "CC":
        vals = ('kcorrRE', 'kcorrIM', 'corrRE', 'corrIM')
        keys = ["/corr/CC_{}*".format(val) for val in vals]

    elif op =="NN":
        vals = ('NN')
        keys = ["/corr/NN_corr*"]

    NNs = [[sorted( glob.glob( file + key), key= lambda x: float(x[ len(file) + len(key)  - 1: ])) for file in files] for key in keys]

    minval = 1e10
    for i, NNkey in enumerate(NNs):
        for j, NNcase in enumerate(NNkey):
            minval = min(minval, len(NNcase))
            
    #print(NNs)
    Corrs = [[[[] for _ in range(minval)] for _ in range(len(NNs[0]))] for _ in range(len(NNs))]


    for i, NNkey in enumerate(NNs):
        for j, NNcase in enumerate(NNkey):
            for k, f in enumerate(NNcase[:minval]):
                raw = np.loadtxt(f)
                if 'kcorr' not in f:
                    raw = raw[2 :-2, 2:-2]

                Corrs[i][j][k] = raw
                

    Corrs = np.array(Corrs)
    lo = np.amin(Corrs, axis=tuple(range(2, Corrs.ndim)))
    hi = np.amax(Corrs, axis=tuple(range(2, Corrs.ndim)))

    print(Corrs.shape, lo.shape, hi.shape)

    fig, axes = plt.subplots(4, Corrs.shape[1], figsize = (5* Corrs.shape[1], 20))

    if Corrs.shape[1] == 1:
        axes = [[ax] for ax in axes]

    axes : list[list[plt.Axes]]
    dividers = [[make_axes_locatable(ax) for ax in row] for row in axes]
    caxes = [[divider.append_axes('bottom', size='10%',  pad=0.05) for divider in row] for row in dividers]

    anim = FuncAnimation(fig, animate, frames=Corrs.shape[2])

    writervideo = animation.FFMpegWriter(fps=10)
    anim.save( outdir + '{}Corr{}.mp4'.format( op, '_'.join(files)), writer=writervideo)

def plot_qe():

    raw = np.loadtxt('12_allenergy')
    ex = raw.shape[-1]
    single_ref = np.arange(0.01, 4.02, 0.01)
    ref = np.repeat(single_ref, ex)
    chain = np.loadtxt('12_chain_energy')
    chain_gap = chain - chain[0]
    #print(raw.shape, ref.shape)
    fig, ax = plt.subplots(figsize = (10, 10))

    ax.set_ylim( np.amin(raw), np.amax(raw))
    ax.scatter(ref, raw.flatten(), s= 0.1)
    #ax.scatter( chain, np.ones( chain.shape[0]) * -19)
    for ene in chain:
        ax.axhline( ene, c='r', linewidth =  0.1)
    ax.scatter( chain_gap, chain, s= 5)
    ax.set_xlabel('QE energy')
    ax.set_ylabel('Excited state energy')
    ax.set_title('Plotting excited states vs. QE energy, 2 QE, L = 12')
    fig.savefig('plots/QEplot12.pdf')

def static_sd(ft=False):

    dirs = sys.argv[2]
    files = glob.glob( dirs + '/occ*')

    times = set()
    Ls = set()

    for f in files:

        timestr, lstr = f.split('_')
        time = float( timestr[ len(dirs) + 1 + len('occ') :])
        L = int( lstr[ len('L') : ])

        times.add(time)
        Ls.add(L)

    fig, axs = plt.subplots(len(Ls), len(times) + 1, figsize = (10 * len(times) , 10 * len(Ls) ))

    if len(files) == 1:
        axs = [[axs]]

    elif len(Ls) == 1:
        axs = [axs]

    
    for r, L in enumerate(sorted(list(Ls))):

        ax : plt.Axes = axs[r][-1]

        freq = np.loadtxt( dirs + '/freq_L{}'.format(L))
        overlap = np.loadtxt( dirs + '/overlap_L{}'.format(L))
        sites = np.loadtxt( dirs + '/sites_L{}'.format(L))

        overlap = np.log10(overlap)
        ax.scatter( freq,overlap)

        for i in range(freq.shape[0]):
            ax.annotate(str(int(sites[i])), xy=( freq[i], overlap[i]))

        ax.set_xlabel('energy')
        ax.set_ylabel(r'$log_{10}$ overlap ')

        for c, t in enumerate(sorted(list(times))):
        
            ax : plt.Axes = axs[r][c]
            exs = np.loadtxt( dirs + '/occ{}_L{}'.format(t, L))
            chains = exs[:, 1:-1]

            ttotal = chains.shape[0]
            ref = np.arange(ttotal) * t
            cmap = colormaps.get_cmap('rainbow', L//2)
            for j in range(L):

                if j < L//2:
                    #marker = '.'
                    linestyle = 'solid'
                else:
                    #marker = 'v'
                    linestyle = 'dotted'

                #ax.scatter( ref, chains[i][:, j], label='site {}'.format(j + 1), marker=marker, s=4)

                ex = chains[:, j]
                if not ft:
                    ax.plot( ref, ex , label='site {}'.format(j + 1), linestyle = linestyle, c=cmap(j%(L // 2)))

                else:
                    X = fft(ex)
                    freq = fftfreq( len(X), t)
                    mid = len(X) // 2 + 1
                    X = X[1:mid]
                    freq = freq[1:mid]

                    pairs = sorted(list(zip(X, freq)), key=lambda x: np.abs(x[0]))[:10]
                    
                    maxfreq = [x[1] for x in pairs]
                    maxval = [np.abs(x[0]) for x in pairs]

                    print(maxfreq)
                    ax.scatter( maxfreq, maxval, label='site {}'.format(j + 1), linestyle = linestyle, c=cmap(j%(L // 2)) )


                handles, labels = ax.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                ax.legend(by_label.values(), by_label.keys())

                if not ft:
                    ax.set_xlim( -10 * t, ttotal * t)
                ax.set_xlabel('time')
                ax.set_ylabel('occ')
                ax.set_title('OCC, time = {}, L = {}'.format(t, L))
                
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    #plt.show()
    fig.savefig( 'plots/sd{}ft{}.pdf'.format(dirs, ft)) 


def occ_direct():

    def animate(i):
        print("frame {}".format(i))
        ax.clear()
        ax.scatter( ref, occ[i], c='blue')

        if systype == "Electron":
            ax.scatter( ref, -occdn[i], c ='red')

        ax.set_ylim( lo, hi)
        ax.set_title('Charge density per site vs. time')

    systype = sys.argv[3]
    file = sys.argv[2]
    outdir = os.getcwd() + '/vids/' 


    if systype == "Electron":
        occup = np.loadtxt(file + '/occup')
        occdn = np.loadtxt(file + '/occdn')
        lo = -np.amax( occdn)    
        hi = np.amax( occup)
        occ = occup
        

    else:
        occ = np.loadtxt(file + '/occ')
        lo = np.amin(occ)
        hi = np.amax(occ)

        
    ref = np.arange(occ.shape[-1])
    
    fig, ax = plt.subplots()
    anim = FuncAnimation(fig, animate, frames=occ.shape[0])

    writervideo = animation.FFMpegWriter(fps=10)
    anim.save( outdir + 'occ_direct{}.mp4'.format( file), writer=writervideo)

def plot_tcd():


    def animate(i):
        print("frame {}".format(i))
        ax.clear()

        for j in range():
            ax.scatter( ref, occup[i], c='blue')
            ax.scatter( ref, -occdn[i], c ='red')

        ax.set_ylim( lo, hi)


    file = sys.argv[2:]
    outdir = os.getcwd() + '/vids/' 

    for n, d in enumerate(file):
        a = 1

    occup = np.loadtxt(file + '/occup')
    occdn = np.loadtxt(file + '/occdn')

    lo = -np.amax( occdn)    
    hi = np.amax( occup)

    ref = np.arange( occup.shape[-1])
    fig, ax = plt.subplots()
    anim = FuncAnimation(fig, animate, frames=occup.shape[0])

    writervideo = animation.FFMpegWriter(fps=10)
    anim.save( outdir + 'occ_direct{}.mp4'.format( file), writer=writervideo)


def plot_gpi():

    dirs = sys.argv[2:]
    
    #cmap =  colormaps.get_cmap('rainbow')

    def animate(i):
        print("frame {}".format(i))

        for s in range(states):
            
            row = s // c
            col = s % c
            ax : plt.Axes = axes[row][col]
            ax.clear()

            #cval = di/plotdim
            ax.scatter( np.arange( sections), gpi[i, s, :])
            ax.plot( np.arange( sections), gpi[i, s, :], label=str(s))
            
            ax.legend()
            ax.set_title('GPI for transition from state {}'.format(s))
            ax.set_xlabel('section')
            ax.set_ylabel('GPI abs')
            ax.set_ylim( lo[s], hi[s])

        fig.tight_layout()


    #fig, axes = plt.subplots(len(dirs), 1, figsize=(30 , 10 * len(dirs)) )
    for _, d in enumerate(dirs):
        
        
        temp_tag = '/gpi_' 
        times  = sorted([ s[ len( d + temp_tag + 'RE') :] for s in glob.glob(d + temp_tag + 'RE*')])

        #print(times)
        states, sections = np.loadtxt( d + temp_tag + 'RE0.5').shape

        r = c = int(np.ceil( np.sqrt(states))) 
        gpi_re = np.zeros((len(times), states, sections))
        gpi_im = np.zeros((len(times), states, sections))

        for i, time in enumerate(times):
            gpi_re[i, :, :] = np.loadtxt( d + temp_tag + 'RE' + time)
            gpi_im[i, :, :] = np.loadtxt( d + temp_tag + 'IM' + time)

        gpi = np.sqrt( np.abs(gpi_re) ** 2  + np.abs(gpi_im) ** 2 )
        
        fig, axes = plt.subplots( r, c, figsize = (5 *c, 4 * r))
        
        lo = [ np.amin(gpi[:, s, :]) for s in range(states)]
        hi = [ np.amax(gpi[:, s, :]) for s in range(states)]

        anim = FuncAnimation(fig, animate, frames=len(times))

        writervideo = animation.FFMpegWriter(fps=10)
        anim.save( 'vids/' + 'gpi_{}.mp4'.format( d), writer=writervideo)


def animate():



    fig, axes  = plt.subplots(2, 1, figsize = (10, 10))
    fig.subplots_adjust(right=0.75)
    L = 100
    xdata = np.arange(L)
    #k, w = 2, 0.001
    #k = [2, 2]
    #w = [0.001, 0.001]
    y = np.zeros((2, L))

    def update_ani(frame):
        
        ax : plt.Axes = axes[0]
        ax.clear()

        y[0] = np.sin( k1.val * xdata + w1.val * frame) * np.exp( - 0.01 * ( xdata -  frame) ** 2) 
        y[1] = np.sin( k2.val * xdata + w2.val * frame) * np.exp( - 0.01 * ( xdata -  frame) ** 2) 

        #y[1] = np.sin( k[1] * (L - xdata) + w[1] * frame) * np.exp( - 0.01 * ( (L- xdata) -  frame) ** 2) 
        ax.plot(xdata, y[0], c='red')
        ax.plot(xdata, y[1], c='blue')
        ax.set_ylim(-1, 1)

        ax : plt.Axes = axes[1]
        ax.clear()

        ax.plot( xdata, y[0] + y[1])
        ax.set_ylim(-1, 1)


    # Make a horizontal slider to control the frequency.
    ax_k1 = fig.add_axes([0.78, 0.1, 0.2, 0.03])
    k1 = Slider(
        ax=ax_k1,
        label='k1',
        valmin=0,
        valmax=30,
        valinit=2
    )

    # Make a vertically oriented slider to control the amplitude
    ax_w1 = fig.add_axes([0.78, 0.2, 0.2, 0.03])
    w1 = Slider(
        ax=ax_w1,
        label="w1",
        valmin=0,
        valmax=1,
        valinit=1
    )

    ax_k2 = fig.add_axes([0.78, 0.3, 0.2, 0.03])
    k2 = Slider(
        ax=ax_k2,
        label="k2",
        valmin=0,
        valmax=30,
        valinit=1
    )

    ax_w2 = fig.add_axes([0.78, 0.4, 0.2, 0.03])
    w2 = Slider(
        ax=ax_w2,
        label="w2",
        valmin=0,
        valmax=1,
        valinit=1
    )


    # The function to be called anytime a slider's value changes
    def update(val):
        
        pass
        #plt.show()
        #line.set_ydata(f(t, amp_slider.val, freq_slider.val))
        #fig.canvas.draw_idle()


    # register the update function with each slider
    k1.on_changed(update)
    w1.on_changed(update)
    k2.on_changed(update)
    w2.on_changed(update)

    ani = FuncAnimation(fig, update_ani, frames= 100,  interval=2)

    plt.show()

def ee():

    files = sys.argv[2:]
    print( [ np.loadtxt('{}/EE'.format(file)).shape for file in files])
    EEs = np.array([ np.loadtxt( '{}/EE'.format(file)) for file in files])
    bonds = np.array([ np.loadtxt( '{}/bonddim'.format(file), dtype=int) for file in files])

    
    SOI = np.loadtxt('SOI')
    fig, axes = plt.subplots( EEs.shape[0], 2, figsize= ( 30, EEs.shape[0] * 10 ))
    cmap = colormaps.get_cmap('rainbow')

    for i, EE in enumerate(EEs):

        points = EE.shape[-1]
        bond = bonds[i]

        ax : plt.Axes = axes[i][0]
        axd : plt.Axes = axes[i][1]
        for j in range(points):

            ax.plot( EE[:, j], label=SOI[j], c = cmap( j /points))
            axd.plot( bond[:, j], label=SOI[j], c = cmap( j /points))

        
        ax.legend()
        ax.set_xlabel('time')
        ax.set_ylabel('Von Neunman Entropy at point')
        ax.set_title(files[i])

        axd.legend()
        axd.set_xlabel('time')
        axd.set_ylabel('Bond dimension at point')
        axd.set_title(files[i])

    fig.savefig('plots/EE.pdf')
        

def current():

    SCALE_STRING = "1"
    SCALE_FACTOR = eval(SCALE_STRING)
    file = sys.argv[2]
    fig, ax = plt.subplots(figsize= (10, 6))
    ax: plt.Axes = ax
    timestep = float(sys.argv[3])

    try:
        raw = np.loadtxt(file + '/current')
        raw = np.insert(raw, 0, 0) 
        ref = np.arange( raw.shape[-1]) * timestep

        ax.plot(ref, raw * SCALE_FACTOR, label= r'Spatial basis: from $ I = -2  \sum_{k} v_k Im \langle  a_k c_S \rangle $')


    except FileNotFoundError:
        print("No current file")
        pass

    try:
        raw = np.loadtxt(file + '/mixcurrent')
        raw = np.insert(raw, 0, 0)
        ref = np.arange( raw.shape[-1]) * timestep

        ax.plot(ref, raw * SCALE_FACTOR, label= r'Mix basis: from $ I = -2  \sum_{k} v_k Im \langle  a_k c_S \rangle $')

    except FileNotFoundError:
        print("No mix current file")
        pass

    try:
        if os.path.exists(file + '/occ'):
            tags = ['']

        else:
            tags = ['up', 'dn']

        for tag in tags:
            rawocc = np.loadtxt(file + '/occ' + tag)
            gs = np.loadtxt(file + '/gs' + tag)

            rawocc = np.insert(rawocc, 0, gs, axis=0)
            LR = np.loadtxt( file + '/LR')

            times, _ = rawocc.shape
            N = LR.shape[0] // 2
            pseudo = np.zeros( times - 1)

            
            ranges = range(N)

            ref = np.arange( pseudo.shape[0]) * timestep
            #calculate the current from occ
            for t in range(1, times - 1):
                pseudo[t ] +=  - (np.sum(rawocc[t + 1, ranges]) - np.sum(rawocc[t - 1, ranges])) / (2 * timestep)

        ax.plot( ref, pseudo * SCALE_FACTOR, label='Spatial basis: From charge density')


    except FileNotFoundError:
        print("No occ file")
        pass

    try:

        if os.path.exists(file + '/mixocc'):
            tags = ['']

        else:
            tags = ['up', 'dn']

        pseudos = []

        for tag in tags:
            rawocc = np.loadtxt(file + '/mixocc' + tag)
            gs = np.loadtxt(file + '/mixgs' + tag)

            rawocc = np.insert(rawocc, 0, gs, axis=0)  
            LR = np.loadtxt( file + '/LR')

            times, alldim = rawocc.shape
            N_LR = LR.shape[0] // 2
            syssize = alldim - N_LR * 2

            print(N_LR, syssize)
            pseudo = np.zeros( times - 1)
            
            
            ids = np.where( LR == 1)[0]

            ranges = [ val + syssize if val >= N_LR else val for val in ids]

            ref = np.arange( pseudo.shape[0]) * timestep
            #calculate the current from occ
            for t in range(1, times - 1):
                pseudo[t ] =  - (np.sum(rawocc[t + 1, ranges]) - np.sum(rawocc[t - 1, ranges])) / ( 2 * timestep)

            pseudos.append( [pseudo])

        ax.plot( ref, np.sum( np.array(pseudos), axis=0).flatten() * SCALE_FACTOR, label='Mix basis: From charge density')


    except FileNotFoundError:
        print("no mix occ file")
        pass

    ax.legend()
    ax.set_title('Current as a function of time, CURRENT_SCALE_FACTOR = ' + (SCALE_STRING))
    ax.set_xlabel('time')
    ax.set_ylabel('I')
    
    fig.tight_layout()
    fig.savefig('plots/current{}.pdf'.format(file))



def time_bond_ee(plot_energy=False, plot_qe=False):

    files = sys.argv[2:]

    if plot_energy:
        fig, axes = plt.subplots( 5, len(files), figsize = (8 * len(files), 20))
        axes = [ [ax] for ax in axes] if len(files) == 1 else axes

    else:
        fig, axes = plt.subplots( 4, len(files), figsize = ( 8 * len(files), 15))
        axes = [[ax] for ax in axes] if len(files) == 1 else axes

    for i, file in enumerate(files):
        
        ee = np.loadtxt(file + '/EE')
        bond = np.loadtxt(file + '/bonds')

        if not plot_qe:
            ee = ee[:, 1:-2]
            bond = bond[:, 1:-2]

        occ = np.loadtxt(file + '/expN')

        
        sites= ee.shape[1]
        time = ee.shape[0]

        print(sites, time)

        if os.path.exists( file + '/timescale' ):
            timescale = np.loadtxt(file + '/timescale')
        
        else:
            timescale = float(file.split('_')[-1])

        if plot_energy:
            energy = np.loadtxt(file + '/mix_basis_energies')

            ax : plt.Axes = axes[-2][i]

            ref = np.arange(1, energy.shape[0] - 1)
            ax.plot(ref, energy[1:-1])
            ax.hlines(0.25, xmin= ref[0], xmax = ref[-1], label='Bias L', colors = 'red')
            ax.hlines(-0.25, xmin= ref[0], xmax = ref[-1], label='Bias R', colors= 'green')
            ax.vlines( energy[ energy < -0.25].shape[0], ymin = energy[0], ymax = energy[-1], colors='orange', linestyles ='dotted')
            ax.vlines( energy[ energy < 0.25].shape[0] + 1, ymin = energy[0], ymax = energy[-1], colors='orange', linestyles ='dotted')
            ax.legend()
        
        ax : plt.Axes =axes[0][i]
        
        xlabel = "Time"
        ylabel = "Site Number"
        extent = [ 1 * timescale, time * timescale, 2, sites + 1]

        titlefile = '_'.join(file.split('_')[:-1]) + ' timescale: {}'.format(timescale) + ' plot_qe: {}'.format(plot_qe)
        im = ax.imshow(bond.transpose(), aspect="auto", extent=extent, interpolation='none')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('Max bond dim. vs. t, each site, {}'.format(titlefile))
        fig.colorbar(im)

        ax : plt.Axes =axes[1][i]

        im = ax.imshow(ee.transpose(), aspect="auto", extent=extent, interpolation='none')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('Von Neumann Entropy on bipartite cut. vs. t. {}'.format(titlefile))
        fig.colorbar(im)

        ax : plt.Axes =axes[2][i]

        seff = np.log( np.power( np.sum( np.exp( 3 * ee) , axis=1) / (ee.shape[-1] -1 ), 1/3))
        times = np.arange( ee.shape[0]) * timescale
        
        ax.plot(times, seff)
        ax.set_xlabel('Time')
        ax.set_ylabel('S_eff')
        ax.set_title('Effectly Entropy vs. t. {}'.format(titlefile))

        ax : plt.Axes =axes[-1][i]

        exl = occ[:, 1]
        exr = occ[:, -2]

        times = np.arange( occ.shape[0]) * timescale
        ax.plot(times, exl, label='QE left')
        ax.plot(times, exr, label='QE right')
        ax.set_xlabel('Time')
        ax.set_ylabel('QE level')
        ax.set_title('QE levels vs. t. {}'.format(titlefile))

        ax.legend()


    fig.tight_layout()
    fig.savefig('plots/{}allEEvtime.pdf'.format('_'.join(files)))



if __name__ == '__main__':
    
    control = int(sys.argv[1])

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

    elif control == 8:
        check_pl_level()

    elif control == 9:
        plot_corr()

    elif control == 10:
        plot_qe()

    elif control == 11:
        static_sd(ft=False)


    elif control == 12:
        plot_gpi()

    elif control == 13:
        occ_direct()

    elif control == 14:
        animate()

    elif control == 15:
        ee()

    elif control == 16:
        current()

    elif control == 17:
        time_bond_ee()