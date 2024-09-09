import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib import rc
from math import ceil
import h5py

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
from visualize import work
import matplotlib as mpl

sys.path.append(os.path.abspath("/Users/knl20/Desktop/Code/TN"))
from dataloader import *




def check_interface(override=''):

    files = sys.argv[1:]
    ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(['effE', 'EE'], files)

    def linear(x, a, b):
        return a * x + b

    fig, axes = plt.subplots(len(files), 1, figsize=(10, 5*len(files)))

    for i, file in enumerate(files):

        offset = int(40/float(override[-3:])) + 5
        ee = ees[i]

        # we cutoff the tail for 0 entanglement


        ee = ee[offset:, :-1]

        logee = np.log(ee)

        diff = logee[:, :-1] - logee[:, 1:]
        diff = np.where(diff > 5, diff, 0)

        diff = np.argwhere(diff)
        ax : plt.Axes = axes[i]

        x = diff[:, 0]
        y = diff[:, 1]

        param, param_cov = curve_fit(linear, x, y)

        a, b = param
        s = 0.4

        ref = np.linspace( np.amin(x), np.amax(x), x.shape[0])
        ax.scatter( x, y, label='extrapolated wavefront', s=s)

        ax.plot( ref, linear(ref, a, b), label='Fit: y = {:.3f}x + {:.3f}'.format(a, b), c='red')
        
        ax.set_ylabel('site number')
        ax.legend()

    fig.savefig( 'plots/{}interface.pdf'.format(override))


def get_qesite(file):

    if 'QEtwo' in file or 'QESSH' in file:
        qeidx = [0, 2]

        total = 2 + 1

    elif 'X' in file:
        qeidx = [
        0 ,
        3 ,
        5 ,
        8 
        ]

        total = 4 + 4 + 1

    else:  
        qeidx = [0, 2, 3, 5]
        total = 2 + 4

    actuals = sorted(list(set(range(total)) - set(qeidx)))

    return qeidx, actuals


def get_plt_attr():
    total = 9

    fig = plt.figure(figsize=(14, 6), layout="constrained")
    spec = fig.add_gridspec(2, 5)

    axes = [[] for _ in range(total)]

    axes[0] = fig.add_subplot(spec[0, 0])
    axes[1] = fig.add_subplot(spec[0, 1])
    axes[2] = fig.add_subplot(spec[0, 3])
    axes[3] = fig.add_subplot(spec[0, 4])
    axes[4] = fig.add_subplot(spec[:, 2])
    axes[5] = fig.add_subplot(spec[1, 0])
    axes[6] = fig.add_subplot(spec[1, 1])
    axes[7] = fig.add_subplot(spec[1, 3])
    axes[8] = fig.add_subplot(spec[1, 4])

    return fig, axes

def get_begin_end(shape, key='parallel'):



    if 'X' in key:

        midpoint = shape//2
        quadpoint = midpoint // 2
        
        begin = {
            0 : 0,
            1 : 2,
            2 : midpoint + 1,
            3 : midpoint + quadpoint - 1,
            4 : midpoint ,
            5 : quadpoint ,
            6 : quadpoint + 2,
            7 : midpoint + quadpoint + 1,
            8 : shape - 2
        }

        end = {
            0 : 2,
            1 : quadpoint,
            2 : midpoint + quadpoint - 1,
            3 : midpoint + quadpoint + 1,
            4 : midpoint + 1,
            5 : quadpoint + 2,
            6 : midpoint,
            7 : shape - 2,
            8 : shape
        }

    elif 'QEtwo' in key:

        #QE1 C QE2
        segment = [0, 2, shape - 2, shape]

        begin = { i : segment[i] for i in range(len(segment) - 1)}
        end = { i : segment[i + 1] for i in range(len(segment) - 1)}

    elif 'DD' in key:

        midpoint = shape//2 - 1

        segment = [0, 2, midpoint - 2, midpoint, midpoint + 2, shape - 4, shape - 2, shape]
        
        begin = { i : segment[i] for i in range(len(segment) - 1)}
        end = { i : segment[i + 1] for i in range(len(segment) - 1)}

    else:
 #QE1, C1, QE2, QE3, C2, QE4

        midpoint = shape//2
        quadpoint = midpoint // 2

        segment = [0, 2, midpoint - 2, midpoint, midpoint + 2, shape - 2, shape]
        
        begin = { i : segment[i] for i in range(len(segment) - 1)}
        end = { i : segment[i + 1] for i in range(len(segment) - 1)}

    return begin, end

def seff(filename_override=''):

    def set_legend(ax : plt.Axes):
        # Shrink current axis by 20%
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.6))
        #ax2.legend(loc='lower right', handles=twins)


    files = sys.argv[1:]

    #seff_override = len(set([ searchkey('U', f) for f in files])) * len(set([ searchkey('tswitch', f) for f in files])) 

    #dimtotal = sorted([ int(searchkey('TEdim', f)) for f in files])

    figseff, axseffs = plt.subplots( figsize=(20, 3))

    ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(['effE'], files)
        

    for i, file in enumerate(files):

        cmap = mpl.cm.hot
        tag = file
        #tag, linestyle= get_tag(file)
        #tag = tag +  'bd = ' + ''.join([ f[len("TEdim"):] if 'TEdim' in f else '' for f in file.split('_')])

        print(tag)
        time = times[i]
        seff : np.ndarray = seffs[i]

        axseff :plt.Axes = axseffs
        color = cmap((i + 1)/( len(files) + 1))
        axseff.plot(time, seff, label='Seff ' + tag, c = color,
                    )
        axseff.set_xlim(tmin, tmax)
        axseff.set_ylim( 0, s_hi)
        axseff.set_xlabel('Time')
        axseff.set_ylabel(r"$S_{eff}$")
        axseff.set_title('Effectly Entropy vs. t')
        #axseff.legend()

    set_legend(axseff)

    figseff.tight_layout()

    if filename_override == '':

        try:
            figseff.savefig( 'plots/' + 'Seff' + '_'.join(files) + '.pdf', dpi=1000)

        except OSError:
            figseff.savefig('plots/test256.pdf')

    else:
        figseff.savefig( 'plots/' + 'Seff' + filename_override + '.pdf', dpi=1000)
    #cnt += 1



def occ_direct(get_animate=True,  file=None, occ=None, fig=None, ax=None, out=True, diff=False):

    def animate(i):
        print("frame {}".format(i))

        ax : plt.Axes = ax

        ref = np.arange( occ.shape[-1])
        ax.clear()
        ax.scatter( ref, occ[i], c='blue')

        # if systype == "Electron":
        #     ax.scatter( ref, -occdn[i], c ='red')

        ax.set_ylim( lo, hi)
        ax.set_title('Charge density per site vs. time')

    #systype = sys.argv[3]

    if not file:
        file = sys.argv[1]

    occ = np.loadtxt(file + '/occ')
    outdir = os.getcwd() + '/vids/' 
    qeidx, actual = get_qesite(file)
    lo = np.amin(occ)
    hi = np.amax(occ)

    begin, end = get_begin_end(occ.shape[-1], key=file)
    ranges = np.concatenate([ np.arange(begin[b], end[b]) for b in actual])

    print(ranges)
    occ = occ[:, ranges]

    if diff:
        occ = occ - occ[0, :]
    
    #fig, axes = get_plt_attr()

    if not fig and not ax:
        fig, ax = plt.subplots(  figsize=(15, 5))

    if get_animate:
        anim = FuncAnimation(fig, animate, frames=occ.shape[0])
        writervideo = animation.FFMpegWriter(fps=10)
        anim.save( outdir + 'occ_direct{}.mp4'.format( file), writer=writervideo)

    else:

        #print(ranges)
        ax : plt.Axes = ax
        im = ax.imshow(occ.transpose(), aspect="auto", cmap='hot',
                       interpolation='none', rasterized=True, 
                       #norm='log',
                       #vmin = 10 ** (-2.5), vmax = 1
                       )
        fig.colorbar(im, location='bottom')

        ax.set_xlabel('time')
        ax.set_ylabel('sites')
        ax.invert_yaxis()

        ax.set_title('Occ vs. time, diff={}'.format(diff))

        if out:
            fig.savefig('plots/occ{}.pdf'.format(file))




def tcd_direct():

    def animate(i):
        print("frame {}".format(i))

        #ax : plt.Axes = ax

        ax : plt.Axes = axes[0]

        ref = np.arange( tcd.shape[-1])
        ax.clear()
        ax.scatter( ref, tcd[i], c='blue')
        ax.plot( ref, tcd[i], c='blue')

        #sep =  np.arange(0, ref.shape[-1], 12)
        #ax.vlines( sep, lo *np.ones(sep.shape[0]), hi *np.ones(sep.shape[0]), linestyles='dotted')

        # if systype == "Electron":
        #     ax.scatter( ref, -occdn[i], c ='red')

        ax.set_ylim( lo, hi)
        ax.set_xlim( 0, tcd.shape[-1])
        ax.set_title('Transition Charge density (TCD) per site vs. time, step = {}'.format(i))

        ax : plt.Axes = axes[1]

        ax.clear()

        ref_processed = np.arange(tcd_processed.shape[-1])
        ax.scatter( ref_processed, tcd_processed[i], c='orange')
        ax.plot( ref_processed, tcd_processed[i], c='orange')
        #ax.vlines( sep, loglo *np.ones(sep.shape[0]), loghi *np.ones(sep.shape[0]), linestyles='dotted')
        # if systype == "Electron":
        #     ax.scatter( ref, -occdn[i], c ='red')

        ax.set_ylim( processed_lo, processed_hi)
        ax.set_xlim( 0, tcd.shape[-1])
        ax.set_title('Transition Charge density (TCD) per site vs. time, LOG, step = {}'.format(i))

    file = sys.argv[1]

    tcd = np.loadtxt(file + '/TCDdyna')

    lo = np.amin(tcd)
    hi = np.amax(tcd)

    # logtcd = np.where( tcd !=0, tcd, np.exp(1e-10))
    # logtcd = np.where( logtcd > 0, -np.log(logtcd), np.log(-logtcd))

    # loglo = np.amin(logtcd)
    # loghi = np.amax(logtcd)

    pad = np.zeros( (tcd.shape[0], 1))
    tcd_processed = np.concatenate( (pad, tcd, pad), axis=1)

    tcd_processed = (tcd_processed[:, 1:] + tcd_processed[:, :-1])/2
    processed_lo = np.amin(tcd_processed)
    processed_hi = np.amax(tcd_processed)

    fig, axes= plt.subplots( 2, 1, figsize=(15, 10))

    outdir = os.getcwd() + '/vids/' 

    anim = FuncAnimation(fig, animate, frames=tcd.shape[0])
    writervideo = animation.FFMpegWriter(fps=2)
    anim.save( outdir + 'tcd_direct{}.mp4'.format( file.split('/')[-1]), writer=writervideo)



def all(override='', tilltick=1000, visualize=True, multisvn =False, diff=False):
    
    def get_seff(ee : np.ndarray):
        return np.log( np.power( np.sum( np.exp( 3 * ee) , axis=1) / (ee.shape[-1] -1 ), 1/3))

    files = sys.argv[1:]

    outdir = os.getcwd() + '/plots/'

    col = 9 - 3 * (not visualize) - 2 * (not multisvn) - (not diff)
    
    fig, axes = plt.subplots(len(files), col, figsize=(15 * col, 5 * len(files)))

    if len(files) == 1:
        axes = [axes]

    ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(['effE', 'EE'], files)


    for i, file in enumerate(files):

        print(file)

        if 'dd' in file:
            key = 'DD'

        elif 'parallel' in file:
            key = 'parallel'

        elif 'QEtwo' in file or 'QESSH' in file:
            key = 'QEtwo'

        else:
            raise(ValueError("Unknown file, recheck"))

        SvN = ees[i]
        occ = occs[i]

        time = times[i]

        print(key)
        print(occ.shape)
        begin, end = get_begin_end(occ.shape[-1], key=key)

        cnt = 0
        ax : plt.Axes = axes[i][cnt]

        occ_direct(get_animate=False, file=file, occ=occ, fig=fig, ax=ax, out=False, diff=False)

        if diff:
            cnt += 1
            ax : plt.Axes = axes[i][cnt]
            occ_direct(get_animate=False, file=file, occ=occ, fig=fig, ax=ax, out=False, diff=True)

        cnt += 1
        ax : plt.Axes = axes[i][cnt]
        xlabel = "Time"
        ylabel = "Site Number"

        actual = file.split('/')[-1]

        # LRBT
        extent = [ 0, time[-1], 1, SvN.shape[-1]]

        im = ax.imshow(SvN.transpose(), aspect="auto", cmap='hot',
                       interpolation='none', rasterized=True, 
                       #norm='log',
                    vmin = ee_lo, vmax=ee_hi, 
                    extent=extent,
                    origin='lower'
                    )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('Von Neumann Entropy on bipartite cut. vs. t. \n {}'.format(actual))
        fig.colorbar(im,  location='bottom')

        cnt += 1

        if multisvn:

            Renyi2 = np.loadtxt(file + '/SRenyi2')
            Renyi3 = np.loadtxt(file + '/SRenyi3')

            ax : plt.Axes = axes[i][cnt]
            xlabel = "Time"
            ylabel = "Site Number"

            im = ax.imshow(Renyi2.transpose(), aspect="auto",  cmap='hot',
                        interpolation='none', rasterized=True,
                        #vmin = ee_lo, vmax=ee_hi, 
                        extent=extent,
                        origin='lower'
                        )
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title('$a=2$ Renyi entropy on bipartite cut. vs. t. \n {}'.format(actual))
            fig.colorbar(im,  location='bottom')

            cnt += 1
            ax : plt.Axes = axes[i][cnt]
            xlabel = "Time"
            ylabel = "Site Number"

            im = ax.imshow(Renyi3.transpose(), aspect="auto", cmap='hot',
                        interpolation='none', rasterized=True,
                        #vmin = ee_lo, vmax=ee_hi, 
                        extent=extent,
                        origin='lower'
                        )
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title('$a=3$ Renyi entropy on bipartite cut. vs. t. \n {}'.format(actual))
            fig.colorbar(im,  location='bottom')

            cnt += 1

        ax : plt.Axes = axes[i][cnt]

        seg, _ = get_qesite(file)

        for j, idx in enumerate(seg):
            
            print(begin[idx])
            qe = occ[:, begin[idx]]
            ax.plot( time, qe, label='QE' + str(j + 1))
            ax.set_xlim(tmin, tmax)

        if key == 'DD':

            dd = occ[:, -2]
            ax.plot( time, dd, label='DD contact' )

        ax2 : plt.Axes = ax.twinx()
        seff = get_seff(SvN)
        ax2.scatter(time, seff, label='Total effective Von Neumann entropy', 
                    marker='x', color='black'
                    #facecolors='none', edgecolors='black'
                    )

        ax.set_title('QE levels and entropy vs. t')
        ax2.set_ylabel('$S_{eff}$')
        ax2.set_ylim(0, s_hi)
        ax2.set_xlim(tmin, tmax)
        ax2.legend()
        ax.legend()

        # cnt += 1

        # ax : plt.Axes = axes[i][cnt]

        # ees = [SvN, Renyi2, Renyi3]
        # labels = ['SvN', 'Renyi2', 'Renyi3']

        # for j, ee in enumerate(ees):
            

        #     ax.plot(time, get_seff(ee), label=labels[j])

        # ax.set_ylim(0, s_hi)
        # ax.set_xlabel('Time')
        # ax.set_ylabel('Total effective entropy')
        # ax.set_title('Total effective entropy in system vs. t')

        if visualize:
            cnt += 1    


            args = {
                "files" : [file + '/HQEdyna'],
                "axes" : [axes[i][cnt:]],
                "fig" : fig
            }

            work(system ="QE", args=args)
        
        # ax : plt.Axes = axes[i][cnt]
        # for j, idx in enumerate([1, 4 ]):
            
        #     mid = (end[idx] - begin[idx]) // 2
        #     seg1 = np.sum(occ[:, begin[idx] : mid] - occ[0, begin[idx] : mid], axis = 1)
        #     ax.plot( time, seg1, label='sub_chain_segment = {}, chain = {}'.format(j * 2 + 1, j + 1 ) )

        #     seg2 = np.sum(occ[:,  mid : end[idx]] - occ[0, mid : end[idx]], axis = 1)
        #     ax.plot( time, seg2, label='sub_chain_segment = {}, chain = {}'.format( (j + 1) * 2, j + 1 ) )


        # ax.set_title("Sum of charge density fluctuation over segments (half of each subchain) ")
        
        ax.legend()

    fig.tight_layout()

    dpi = 300
    if override:
        fig.savefig(outdir + '{}.pdf'.format(override), dpi=dpi)
    else:
        fig.savefig(outdir + 'EE{}.pdf'.format('_'.join(files)), dpi=dpi)


def corr():

    filename = "/Users/knl20/Desktop/Code/TN/src/work/corr.h5"

    with h5py.File(filename, "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        print("Keys: %s" % f.keys())


def segment(override=''):

    files = sys.argv[1:]


    segs = [10, 20, 30, 34, 35, 36, 40]
    
    fig, axes = plt.subplots(len(files), len(segs) + 4, figsize=(5 * (len(segs) + 4), 4 * len(files)))
    minval = 10
    maxval = 0

    
    for i, file in enumerate(files):

        args = {
            "files" : [file + '/HQEdyna'],
            "axes" : [axes[i][len(segs) + 1:]],
            "fig" : fig
        }

        work(system ="QE", args=args)

        occ = np.loadtxt(file + '/occ')
        times = np.loadtxt(file + '/times')

        if 'dd' in file:
            key = 'DD'

        elif 'parallel' in file:
            key = 'parallel'

        elif 'QEtwo' in file:
            key = 'QEtwo'

        else:
            key = 'parallel'

        begin, _ = get_begin_end(occ.shape[-1], key=key)
        actual = file.split('/')[-1]

        for j, seg in enumerate(segs):

            ax :plt.Axes = axes[i][j]

            upstart, upend = begin[1], begin[1] + seg
            

            UPLEFT = np.sum( occ[:, upstart : upend], axis=1)
            DUP = (UPLEFT[2:] - UPLEFT[:-2]) /2
            minval = min(minval, np.amin(DUP))
            maxval = max(maxval, np.amax(DUP))
            ax.plot( times[1:-1],  DUP, label='upper {} to {}'.format(upstart + 1, upend ))

            if 'QEtwo' not in file:

                lostart, loend = begin[4], begin[4] + seg
                LOLEFT = np.sum( occ[:, lostart : loend], axis=1)
                DLO = (LOLEFT[2:] - LOLEFT[:-2]) /2
                minval = min(minval, np.amin(DLO))
                maxval = max(maxval, np.amax(DLO))
                ax.plot( times[1:-1], DLO, label='lower {} to {}'.format(lostart + 1, loend ))

            ax.set_title(actual + ', sum 1 to ' + str(seg))
            ax.set_xlabel('time')
            ax.set_ylabel(r'$ \langle dN_L/dt\rangle$')
            ax.legend()

        ax : plt.Axes = axes[i][ len(segs) ]

        qeidx, _ = get_qesite(file)

        for j, idx in enumerate(qeidx):
            
            qe = occ[:, begin[idx]]
            ax.plot( times, qe, label='QE' + str(j + 1))

        ax.set_title('QE levels')
        ax.legend()

    for i in range(len(files)):
        for j in range(len(segs)):

            axes[i][j].set_ylim(minval, maxval)
        
        
    fig.tight_layout()
    fig.savefig('plots/DNDt{}.pdf'.format(override))


if __name__ == '__main__':

    
    override= 'GaussianQE100dyna'
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 
    #occ_direct(get_animate=False)

    #tcd_direct()
    #check_interface(override=override)
    # seff(filename_override=override)
    all(override=override, 
        visualize=False,
        multisvn=False
        #tilltick=150
        )
    #segment(override = override )
    #corr()
