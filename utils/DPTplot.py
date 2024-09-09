import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import h5py
from matplotlib import rc
import matplotlib as mpl
from math import ceil
import time
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
from searchkey import *
from strsort import *
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial
from itertools import product


sys.path.append(os.path.abspath("/Users/knl20/Desktop/Code/TN"))
from dataloader import *

import h5py


def MFloader(L, U, tswitch):

    mfstr = '/Users/knl20/Desktop/Code/TN/LT/meanfield/Current_dpt_mf_N{}_mu0.5_U{}_tswitch{}'.format(L, U, tswitch)
    try:
        raw = np.loadtxt(mfstr)
    except FileNotFoundError:
        srchstr = mfstr +'*'

        print(srchstr)
        f = glob.glob(srchstr)[0]
        raw = np.loadtxt(f)

    return raw

def str_generator(files, colorgroup, linestylegroup, markergroup, exclude_tag={}):

    #tags = sorted([ get_tag(file) for file in files])

    tags = [ get_tag(file, exclude_tag=exclude_tag) for file in files]

    cval_ref = np.linspace(0.1, 0.9, len(tags))
    linestyle_str = [
     'solid',      # Same as (0, ()) or '-'
     'solid',    # Same as (0, (1, 1)) or ':'
     'dashed',    # Same as '--'
     'dashdot']  # Same as '-.'
    
    marker_str = ['o', 'x', '^']

    comb = product(linestyle_str, marker_str)

    if colorgroup == {}:
        cvals = cval_ref

    else:

        cvals = []

        for file in files:
            
            cval = cval_ref[0]
            for key in colorgroup:
                if key in file:
                    cval = colorgroup[key]

            cvals += [cval]


    if linestylegroup == {}:
        linestyles = [ val[0] for val in comb]

    else:

        linestyles = []

        for file in files:
            
            linestyle = 'solid'
            for key in linestylegroup:
                if key in file:
                    linestyle = linestylegroup[key]

            linestyles += [linestyle]


    if markergroup == {}:
        markers = [ val[1] for val in comb] 

    else:

        markers = []

        for file in files:

            marker = ''
            for key in markergroup:
                if key in file:
                    marker = markergroup[key]

            markers += [marker]


    return tags, cvals, linestyles, markers


def get_tag(file, exclude_tag={'ddposition', 'avg'}):

    tag = {}

    if 'basis' not in exclude_tag:
        if 'mix' not in file or 'mixedFalse' in file:
            tag['basis'] = 'spatial'
            #tag['includeQPC'] = 'T'
            

        else: 

            if 'random' in file:
                tag['basis'] = 'mixed random'

            elif 'orderingABSSORTED' in file:
                tag['basis'] = 'mixed abs ordered'

            #elif 'orderingSORTED' in file:
            else:
                tag['basis'] = 'mixed ordered'

            if 'includeUFalse' in file:
                tag['includeQPC'] = 'excludeQPC'

            else:
                tag['includeQPC'] = 'includeQPC'

    if 'ddposition' not in exclude_tag:
        if 'ddpositionL' in file:
            tag['DDpos'] = 'L'

        elif 'ddpositionR' in file:
            tag['DDpos'] = 'R'

        elif 'ddpositionM' in file:
            tag['DDpos'] = 'M'

    if 'avg' not in exclude_tag:
        if 'avgTrue' in file:
            tag['avg'] = 'no DD'

        else:
            tag['avg'] = 'DD'

    if 'dim' not in exclude_tag:
        tag['dim'] = '$\dim = {}$'.format(searchkey( 'TEdim', file))

    if 'L' not in exclude_tag:
        tag['L'] = '$L = {}$'.format(searchkey('L', file))
    #tag['class'] = 'best'
    #tag['L'] = 'dim ' + searchkey('dim', file)
    #res = ', '.join([ key + '=' + tag[key] for key in tag])
    res = ', '.join( tag.values())

    
    return res




def occ_direct():

    def animate(i):
        print("frame {}".format(i))

        ax.clear()
        ax.scatter( ref, occ[i], c='blue')

        # if systype == "Electron":
        #     ax.scatter( ref, -occdn[i], c ='red')

        ax.set_ylim( lo, hi)
        ax.set_title('Charge density per site vs. time')

    #systype = sys.argv[3]
    file = sys.argv[1]
    outdir = os.getcwd() + '/vids/' 


    # if systype == "Electron":
    #     occup = np.loadtxt(file + '/occup')
    #     occdn = np.loadtxt(file + '/occdn')
    #     lo = -np.amax( occdn)    
    #     hi = np.amax( occup)
    #     occ = occup
        

    # else:
    occ = np.loadtxt(file + '/occ')
    lo = np.amin(occ)
    hi = np.amax(occ)

        
    ref = np.arange(occ.shape[-1])
    
    fig, ax = plt.subplots()
    anim = FuncAnimation(fig, animate, frames=occ.shape[0])

    writervideo = animation.FFMpegWriter(fps=10)
    anim.save( outdir + 'occ_direct{}.mp4'.format( 'test'), writer=writervideo)


def corr_direct(override='', get_animate=False, twoD=True, inputoverride=[]):

    def animate(j, axis_offset =0, axis_total=1):

        print("frame {}".format(j))

        for i, file in enumerate(files):

            ax : plt.Axes = axes[i][axis_offset]
            cax : plt.Axes = caxes[i][axis_offset]

            ax.clear()
            cax.clear()
            im = ax.imshow( np.real(corrs[i][keys[i][j]]), cmap='hot_r', vmin=max(lo_cc[i], 1e-8), vmax=hi_cc[i], norm='log')
            ax.invert_yaxis()

            
            if float(keys[i][j]) < 0:
                stage = 1

            elif float(keys[i][j]) < tswitch[i]:
                stage = 2

            else:
                stage = 3

            ax.set_title(r'$|\langle c_i c_j \rangle|$, ' +  '{}, t = {}, Stage {}'.format(tags[i], keys[i][j], stage))

            fig.colorbar(im, cax=cax)

            ax : plt.Axes = axes[i][axis_offset + axis_total]
            cax : plt.Axes = caxes[i][axis_offset + axis_total]

            ax.clear()
            cax.clear()

            nn = np.abs(nns[i][keys[i][j]])

            occ = occs[i][j]
            mask = np.outer(occ, occ)

            truenn = nn - mask
            im = ax.imshow( truenn, cmap='hot_r', vmin=max(lo_nn[i], 1e-8), vmax=hi_nn[i], norm='log')         
            ax.invert_yaxis()

            ax.set_title(r'$\langle n_i n_j \rangle - \langle n_i\rangle \langle n_j\rangle $, ' + '{}, t = {}, Stage {}'.format(tags[i], keys[i][j], stage))

            fig.colorbar(im, cax=cax)


        fig.tight_layout()

    if inputoverride:
        files = inputoverride

    else:
        files = sys.argv[1:]
    outdir = os.getcwd() + '/vids/' 

    tswitch = [float(searchkey('tswitch', file)) for file in files]


    occs = [ np.loadtxt(file + '/occ') for file in files]
    tags = [ get_tag(file)[0] + ', U = ' + searchkey('U', file)[-3:] + ', dim = ' + searchkey('TEdim', file) for file in files]

    corrs = [h5py.File(file + '/corrCdagC.h5','r') for file in files]
    nns = [h5py.File(file + '/corrNN.h5','r') for file in files]
    frames = np.amin( [ len(corr.keys()) for corr in corrs])
    
    keys = [sorted( corr.keys(), key=lambda x: float(x)) for corr in corrs]
    
    # lo_re = [np.amin( np.real([corr[key] for key in corr.keys()])) for corr in corrs]
    # lo_im = [np.amin( np.imag([corr[key] for key in corr.keys()])) for corr in corrs]

    lo_cc = [np.amin( np.abs([corr[key] for key in corr.keys()])) for corr in corrs]
    lo_nn = [np.amin( [nn[key] for key in nn.keys()]) for nn in nns]

    # hi_re = [np.amax( np.real([corr[key] for key in corr.keys()])) for corr in corrs]
    # hi_im = [np.amax( np.imag([corr[key] for key in corr.keys()])) for corr in corrs]

    hi_cc = [np.amax( np.abs([corr[key] for key in corr.keys()])) for corr in corrs]
    hi_nn = [np.amax( [nn[key] for key in nn.keys()]) for nn in nns]

    lo_nn = np.real(lo_nn)
    hi_nn = np.real(hi_nn)
    

    if get_animate:

        fig, axes = plt.subplots(len(files), 2, figsize=(30, 10 * len(files)))

        if len(files) == 1:
            axes = [axes]

        dividers = [[make_axes_locatable(ax) for ax in row] for row in axes]
        caxes = [[divider.append_axes('right', size='10%',  pad=0.05) for divider in row] for row in dividers]
        anim = FuncAnimation(fig, animate, frames=frames)
        writervideo = animation.FFMpegWriter(fps=5)
        anim.save( outdir + 'corr_direct{}.mp4'.format(override ), writer=writervideo)
        #writervideo = animation.HTMLWriter(fps=10)
        #anim.save( outdir + 'corr_direct{}.html'.format(override ), writer=writervideo)

    elif not twoD:

        times = ['0.0', '10.0', '20.0', '30.0']
        fig, axes = plt.subplots(len(files), len(times) * 2, figsize=(20 * len(times), 6 * len(files)))

        if len(files) == 1:
            axes = [axes]

        for i, file in enumerate(files):

            for j, time in enumerate(times):

                s = 2
                CUTOFF = 200
                
                corr = np.array(corrs[i][time]).flatten()
                re = sorted(np.real(corr))[::-1]
                im = sorted(np.imag(corr))[::-1]
                nn = sorted(np.real(np.array( nns[i][time]).flatten()))[::-1]

                ax:plt.Axes = axes[i][j]
                ref = np.arange(1, CUTOFF + 1)
                ax.vlines([35], ymin=-0.3, ymax=1.1, linestyles='dashed', color='red' )
                ax.scatter( ref, re[:CUTOFF], label='real', s=s)
                ax.scatter( ref, im[:CUTOFF], label='imag', s=s)
                #ax.scatter( ref, nn[:CUTOFF], label='NN', s=s)
                ax.set_title('First ' + str(CUTOFF) + ', T = ' + time + ' ' +  tags[i] + ' U = ' + searchkey('U', file)[-3:])
                ax.legend()

                ax:plt.Axes = axes[i][j + len(times)]
                ref = np.arange( len(re) - CUTOFF + 1, len(re) + 1 )
                ax.scatter( ref, re[-CUTOFF:], label='real', s=s)
                ax.scatter( ref, im[-CUTOFF:], label='imag', s=s)
                #ax.scatter( ref, nn[-CUTOFF:], label='NN', s=s)
                ax.set_title('Last ' + str(CUTOFF) + ', T = ' + time + ' ' +  tags[i] + ' U = ' + searchkey('U', file)[-3:])
                ax.legend()
        fig.savefig( 'plots/{}corrstatic.pdf'.format(override))


    else:


        times = [ '0.0', '8.0'] + [ str(s) for s in np.arange(16.0, 36.0, 4.0)]

        fig, axes = plt.subplots(len(files), 2 * len(times), figsize=(2 * 8 * len(times), 10 * len(files)))

        if len(files) == 1:
            axes = [axes]

        dividers = [[make_axes_locatable(ax) for ax in row] for row in axes]
        caxes = [[divider.append_axes('right', size='10%',  pad=0.05) for divider in row] for row in dividers]

        for k, t in enumerate(times):

            idx = keys[0].index(t)

            animate(idx, axis_offset=k, axis_total=len(times))

        fig.savefig( 'plots/{}2Dcorrstatic.pdf'.format(override))




def current_diff(inputoverride=[]):

    def set_color(cur_dim):

        alldim = sorted(list(set([ int(searchkey('TEdim', f)) for f in files])))
        idx = alldim.index( int(cur_dim))

        return idx/len(alldim)
        
    if inputoverride:
        files = inputoverride

    else:
        files = sys.argv[1:]

    ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(['current'], files)

    maxdim = str(max([ int(searchkey('TEdim', f)) for f in files]))

    Us = set([searchkey('U', file) for file in files])

    fig, axes = plt.subplots(len(Us), 2, figsize=(20, 6))

    axdict = { U : i for i, U in enumerate(Us)}

    if len(Us) == 1:
        axes = [axes]

    cmap = mpl.cm.hot

    for i, file in enumerate(files):
        tag = get_tag(file)

        U = searchkey('U', file)

        tag = tag +  'bd = ' + ''.join([ searchkey('TEdim', f) for f in file.split('_')])
    
        curdim = searchkey('TEdim', file)
        replaced = file.replace(curdim, maxdim)
        ref = np.loadtxt( replaced + '/currentLR')

        # if 'abs' not in tag and '32' not in file:
        #     ref = np.concatenate( ([0.0 for _ in range(12)], ref))

        time = times[i][:ref.shape[0]]

        cval = set_color(curdim)
        

        if replaced != file:
            current = currents[i]

            ax : plt.Axes = axes[axdict[U]][0]

            div = np.where( ref > 1e-6, ref, 1)
            diff = (current[:ref.shape[0]] - ref)/div
            ax.plot( time, diff, label=tag, c=cmap(cval), linestyle='dotted')

            ax: plt.Axes = axes[axdict[U]][1]

            diff = np.abs(diff)
            ax.plot( time, diff, label=tag, c=cmap(cval), linestyle='dotted')

    tswitch = searchkey('tswitch', files[0])

    for U, i in axdict.items():

        ax : plt.Axes = axes[i][0]
        ymin, ymax = ax.get_ylim()
        ax.vlines([0, tswitch], ymin, ymax, linestyles='dotted')
        ax.set_xlabel('Time')
        ax.set_title('$L_2$ current difference, reference to $dim = 1024, U = $' + U)
        ax.set_ylabel('$ (I - I_{1024})/I_{1024}$')
        ax.legend()

        ax : plt.Axes = axes[i][1]
        ymin, ymax = ax.get_ylim()
        ax.vlines([0, tswitch], ymin, ymax, linestyles='dotted')
        ax.set_xlabel('Time')
        ax.set_title('$L_2$ current difference, reference to $dim = 1024, U = $' + U)
        ax.set_ylabel('$|(I - I_{1024})/I_{1024}|$')
        ax.set_yscale('log')
        ax.legend()
    #plt.show()

    fig.savefig('plots/currentcomptest.pdf')
    #return fig




# def seff(inputoverride=[], filename_override=''):

#     def set_legend(ax : plt.Axes):
#         # Shrink current axis by 20%
#         ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.6))
#         #ax2.legend(loc='lower right', handles=twins)

#     if inputoverride:
#         files = inputoverride

#     else:
#         files = sys.argv[1:]


#     seff_override = len(set([ searchkey('U', f) for f in files])) * len(set([ searchkey('tswitch', f) for f in files])) * len(set([ searchkey('L', f) for f in files]))

#     figseff, axseffs = plt.subplots(seff_override, 3, figsize=(30, 5 * seff_override))

#     if seff_override == 1:
#         axseffs = [axseffs]

#     ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(['effE', 'EE'], files)
            
#     each = len(files) // seff_override
#     U = None
#     tswitch = None
#     new = False
#     imgcnt = -1

#     for i, file in enumerate(files):

#         cmap = mpl.cm.hot
#         cur_U = searchkey('U', file)
#         cur_tswitch = float(searchkey('tswitch', file))

        
#         tag, linestyle= get_tag(file)
#         tag = tag +  'bd = ' + ''.join([ f[len("TEdim"):] if 'TEdim' in f else '' for f in file.split('_')])

#         print(tag)
#         cur_dim = searchkey('TEdim', file)
#         cur_pos = searchkey('ddposition', file)
                     
#         if cur_U != U or cur_tswitch != tswitch:
            
#             if U or tswitch:
#                 set_legend(axseff)

            
#             tswitch = cur_tswitch
#             U = cur_U
#             new = True
#             internal_cnt = 0

#             imgcnt += 1
        
#         else:
#             new = False
#             internal_cnt += 1

#         verts = [0, tswitch]

#         print(imgcnt, new, internal_cnt)
#         time = times[i]
#         seff : np.ndarray = seffs[i]

#         cnt = 0
#         axseff :plt.Axes = axseffs[imgcnt][cnt]

#         color = cmap(cval)
#         axseff.plot(time, seff, label='Seff ' + 'U={}, '.format(cur_U) + tag, linestyle=linestyle, c = color,
#                     )
#         axseff.set_xlim(tmin, tmax)
#         axseff.vlines(verts, [0] * len(verts), [10] * len(verts), linestyles='dashed')
#         axseff.set_ylim( 0, s_hi)
#         axseff.set_xlabel('Time')
#         axseff.set_ylabel(r"$S_{eff}$")
#         axseff.set_title('Effectly Entropy vs. t, U = ' + str(U)   +', $t_{switch}$ = ' + str(tswitch) )

#         bipart2 = None
#         ee : np.ndarray = ees[i]

        
#         if 'ddpositionL' in file:
#             ddbtw = 0
#             bipart = ddbtw + 1

#         elif 'ddpositionM' in file:
#             ddbtw = ee.shape[-1]//2 
#             bipart = ddbtw - 1
#             bipart2 = ddbtw + 1

#         elif 'avgTrue' in file:
#             ddbtw = ee.shape[-1]//2 -1 
#             bipart = ddbtw - 1
#             bipart2 = ddbtw + 1

#         else:
#             ddbtw = -2
#             bipart = ddbtw - 1

#         print(f'EE shape: {ee.shape}')

#         cnt += 1
#         axseff : plt.Axes = axseffs[imgcnt][cnt]
#         axseff.set_ylim( -0.1, ee_hi)
#         axseff.vlines(verts, [0] * len(verts), [10] * len(verts), linestyles='dashed')
#         axseff.plot(time, ee[:, ddbtw], label= tag + ' b/w two level', linestyle=linestyle, c = color,)
#         axseff.set_title('EE between two levels')

#         if new:
#             m = 'x'
#             s = 0.1
#             axseff.scatter(time, np.ones(time.shape[0]) * np.log(2), label='ln(2)', s=s, marker=m)
#             axseff.scatter(time, np.ones(time.shape[0]) * 0.5 * np.log(2), label='ln(2)/2',  s=s, marker=m)
#             axseff.scatter(time, np.ones(time.shape[0]) * 1.5 * np.log(2), label='3ln(2)/2',  s=s, marker=m)

#         cnt += 1
#         axseff : plt.Axes = axseffs[imgcnt][cnt]

#         if new:
#             m = 'x'
#             s = 0.1
#             axseff.scatter(time, np.ones(time.shape[0]) * np.log(2), label='ln(2)', s=s, marker=m)
#             axseff.scatter(time, np.ones(time.shape[0]) * 0.5 * np.log(2), label='ln(2)/2',  s=s, marker=m)
#             axseff.scatter(time, np.ones(time.shape[0]) * 1.5 * np.log(2), label='3ln(2)/2',  s=s, marker=m)
        
#         axseff.set_ylim( -0.1, ee_hi)
#         axseff.vlines(verts, [0] * len(verts), [10] * len(verts), linestyles='dashed')

#         # if 'ddpositionM' in file or 'avgTrue' in file:
#         #     axseff.plot(time, ee[:, bipart] - ee[:, ddbtw], label= tag + 'b/w DD, sys (adjusted)', linestyle=linestyle, c = color,)
#         # else:
#         axseff.plot(time, ee[:, bipart], label= tag + ' b/w DD and sys', linestyle=linestyle, c = color,)

#         #if bipart2:
#         #    axseff.scatter(time, ee[:, bipart], label= tag + 'btw DD, sys 2', marker='o', c = color,)
#         axseff.set_title('EE on bipartite cut')
#         #axseff.legend()

#     set_legend(axseff)

#     figseff.tight_layout()

#     if filename_override == '':

#         try:
#             figseff.savefig( 'plots/' + 'Seff' + '_'.join(files) + '.pdf', dpi=1000)

#         except OSError:
#             figseff.savefig('plots/test256.pdf')

#     else:
#         figseff.savefig( 'plots/' + 'Seff' + filename_override + '.pdf', dpi=1000)
#     #cnt += 1

#     return figseff




def time_bond_ee( inputfiles : dict = {},  filetag = 'multi', fit=False, save_individual=False, extend_left=False, shadeoff=0, MF=False ):


    def plotcurrent():

        cnt = cat.index('current')
        current = currents[j]
        ax : plt.Axes =axes[ r + cnt][c]
        color = currentcmap(cval)

        if extend_left:
            if time[0] > -tswitch:
                ctime = np.concatenate(( [-tswitch], time))
                current = np.concatenate( ([0], current))

        else:
            ctime = time
        
        if not current_occ_separate:
            axiscolor = color
            ax.tick_params(axis='y', labelcolor=color)

        else:
            axiscolor='black'

        ax.plot(ctime, current, label=tag, color=color, 
            linestyle=linestyle,
            marker=marker, markersize=markersize, linewidth=linewidth,
            )
        
        if len(yscale_current):
            ax.vlines([0, tswitch],  yscale_current[0] * np.ones(2), yscale_current[1] * np.ones(2), linestyles='dotted')

        else:
            ax.vlines([0, tswitch],  np.amin(current) * np.ones(2), np.amax(current) * 1.1 * np.ones(2), linestyles='dotted')



        if cnt == len(cat) - 1:

            ax.set_xticks([-tswitch, 0, tswitch, 2* tswitch])
            ax.set_xlabel(r'$t \ (1/\omega_0)$',  fontsize=fontsize)
        
        else:
            ax.set_xticklabels([])

        if extend_left:
            ax.set_xlim( -tswitch, tswitch*2)

        if not len(yscale_current):
            ax.annotate('* $\sim 1/50$ uniform scale for $U=2.0, 3.0$', xy = (0.02, 0.9), xycoords='axes fraction', fontsize=fontsize/2)


        ax.set_ylabel('Current' + ('*' if not len(yscale_current) else '' ) + ' $(\omega_0)$', color=axiscolor,  fontsize=fontsize)
        ax.set_title(r'$I_{{L\rightarrow R}}$ vs. t {}'.format(group_id), x=xtitle, y=ytitle, fontsize=fontsize)

    def plot_reference_current(L, U, linestyle, marker):

        
        ax : plt.Axes =axes[ r + cat.index('current')][c]
        refcolor='green'
        ref = np.loadtxt('/Users/knl20/Desktop/Code/TN/non-interacting/{}currentCC{}'.format(L, U))

        reftime = np.arange(0, tswitch, 0.25)
        ax.plot( reftime, ref[:reftime.shape[0]], label=r'$I_{L\rightarrow R}$, Stage 2, exact, '+ '$L={}$'.format(L), color=refcolor,
                linestyle=linestyle, markersize=markersize, marker='x', linewidth=linewidth
                )

        if MF:
            refcolor='violet'
            refraw = MFloader(L, U, tswitch)

            reftime = refraw[:, 0]
            ref = refraw[:, 3]

            ax.plot( reftime, ref[:reftime.shape[0]], label=r'$I_{L\rightarrow R}$, MF, ' + '$L={}$'.format(L), color=refcolor,
                    linestyle=linestyle, marker=marker, markersize=markersize, linewidth=linewidth,
                    )

    def plot_reference_occ(L, U, linestyle, marker):

        if not current_occ_separate:
            ax2 : plt.Axes = axes[ r + cat.index('current')][c].twinx()

        else:
            ax2 : plt.Axes = axes[r+ cat.index('occ')][c]

        if MF:
            refcolor='violet'
            refraw = MFloader(L, U, tswitch)

            reftime = refraw[:, 0]
            ref = refraw[:, 5]

            ax2.plot( reftime, ref[:reftime.shape[0]], label='lower dot density, MF, $L={}$'.format(L), color=refcolor,
                    linestyle=linestyle, marker=marker, markersize=markersize, linewidth=linewidth,
                    )



    def plotocc():
        occ = occs[j]
        color = occcmap(cval)

        if not current_occ_separate:
            cnt = cat.index('current')
            ax2 : plt.Axes = axes[ r + cnt ][c].twinx()

        else:
            cnt = cat.index('occ')
            ax2 : plt.Axes = axes[r+ cnt][c]
        

        if 'ddpositionL' in file:
            ddlow = 0

        elif 'ddpositionM' in file or 'avgTrue' in file: 
            ddlow = occ.shape[-1]//2 - 1

        else:
            ddlow = -2

        LT1 = occ[:, ddlow]

        if extend_left:
            if time[0] > -tswitch:
                ctime = np.concatenate(( [-tswitch], time))
                LT1 = np.concatenate( ([1], LT1))

        else:
            ctime = time
            
        obj,  = ax2.plot(ctime, LT1, label=tag, color=color,
        linestyle=linestyle, marker=marker, markersize=markersize, linewidth=linewidth,
        )

        if not current_occ_separate:
            axiscolor = color
            ax2.tick_params(axis='y', labelcolor=color)

        else:
            axiscolor='black'

        
        ax2.set_ylabel(r'$\langle n_1\rangle$', color=axiscolor,  fontsize=fontsize)

        if cnt == len(cat) - 1:
            ax2.set_xlabel('t $(1/\omega_0)$',  fontsize=fontsize)

        else:
            ax2.set_xticklabels([])

        if extend_left:
            ax2.set_xlim( -tswitch, tswitch*2)

        ax2.vlines([0, tswitch], np.zeros(2), 1.1 * np.ones(2), linestyles='dotted')
        ax2.set_title(r'$\langle n_1\rangle$ vs. t {}'.format(group_id), fontsize=fontsize, x=xtitle, y=ytitle)

    def plotee():
        ee = ees[j]
        cnt = cat.index('EE')
        ax : plt.Axes =axes[r + cnt][c]

        extent = [ tmin, tmax, 1, ee.shape[-1]]
        xlabel = "Time"
        ylabel = "Site Number"

        #extent = [ 0, time[-1], 1, ee.shape[-1]]

        im = ax.imshow(ee.transpose(), aspect="auto", extent=extent, cmap='hot',
                    interpolation='none', rasterized=True,
                    #vmin = ee_lo, vmax=ee_hi, 
                    origin='lower'
                    )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('Von Neumann Entropy on bipartite cut. vs. t. \n {}'.format(group_id))
        fig.colorbar(im,  location='bottom')
 
    def plotbond():
        cnt = cat.index('bond')
        ax : plt.Axes =axes[r + cnt][c]
        bond = bonds[j]

        extent = [ tmin, tmax, 1, bond.shape[-1]]
        
        xlabel = "Time"
        ylabel = "Site Number"

        im = ax.imshow(bond.transpose(), aspect="auto", extent=extent, cmap='hot',
                    interpolation='none', rasterized=True,
                    vmin=bd_lo, vmax=bd_hi,
                    origin='lower'
                    )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('Max bond dim. vs. t, each site, \n {}'.format(group_id))
        fig.colorbar(im,  location='bottom')

    def plotseff():
        seff = seffs[j]
        cnt = cat.index('effE')

        ax : plt.Axes =axes[ r + cnt][c]
        color = currentcmap(cval)

        if extend_left:
            if time[0] > -tswitch:
                ctime = np.concatenate(( [-tswitch], time))
                seff = np.concatenate( ([seff[0]], seff))

        else:
            ctime = time


        ax.plot(ctime, seff, label=tag, color=color, 
            linestyle=linestyle,
            marker=marker, markersize=markersize, linewidth=linewidth,
            #marker='x', s=5
            )
        
        ax.vlines([0, tswitch], np.zeros(2), 1.1 * np.amax(seff) * np.ones(2), linestyles='dotted')

        if cnt == len(cat) - 1:
            ax.set_xticks([-tswitch, 0, tswitch, 2* tswitch])
            ax.set_xlabel(r'$t \ (1/\omega_0)$',  fontsize=fontsize)

        else:
            ax.set_xticklabels([])

        if extend_left:
            ax.set_xlim( -tswitch, tswitch*2)

        ax.set_ylabel('Effective Entropy',  fontsize=fontsize)
        ax.set_title(r'$\ln \sqrt{ \frac{ \sum_i^L e^{S_i} }{L} } $' + ' vs. t {}'.format(group_id), fontsize=fontsize, x=xtitle, y=ytitle)

    # input files are grouped by dicts, 
    # each dict: 
    # input = { "page1 title": {"group 1 on page1" : [], "group 2 on page1 : []", etc.}
    #

    def fitocc():

        #def f(x, A, B, w, p, N):
            #return (A + B * np.cos( w * (x - tswitch)) )/ (x - tswitch) ** p + N
            #return n*+[A+B*cos(w*t)]/t^p

        def f(x, A, p, N):
            #return N + A/ (x - tswitch ) **p
            return N + A * np.exp( -p * (x-tswitch))
        
        occ = occs[0]
        LT = occ[:, -2]
        time = times[0]

        if searchkey('U', group[0]) == '3.0':
            #toff = int(searchkey('L', group[0]))/5
            toff = int(searchkey('L', group[0]))/10
        else:
            toff = 0

        idx = np.argwhere( time > (tswitch + toff )).flatten()
        to_fit = LT[idx]
        t_sec = time[idx]

        
        #p0 = (0.1, 1, 1, 2, 0.5)

        if searchkey('U', group[0]) == '5.0':
            p0 = (0.98)

            coeff, pcov = curve_fit( lambda x, a: a, t_sec, to_fit, p0=p0)
            cov = np.sqrt(np.diag(pcov))

        else:

            if searchkey('U', group[0]) == '3.0':
                p0 = (1.8, 0.7, 0.3)

            else:
                p0 = (0.5, 0.1, 0.4)

            coeff, pcov = curve_fit(f, t_sec, to_fit, p0=p0)
            cov = np.sqrt(np.diag(pcov))

        if not current_occ_separate:
            ax2 : plt.Axes = axes[ r + cat.index('current')][c].twinx()

        else:
            ax2 : plt.Axes = axes[r+ cat.index('occ')][c]


        #ax2.plot( t_sec, f(t_sec, *coeff), label=r"FIT STAGE 3 DD : $  \frac{{ {:2f} +  {:2f} * \cos( {:2f} t) }}{{ t^{{{:2f}}} }} + {:2f} $".format(*coeff), linewidth=linewidth, c='black')
        ax2.scatter( tswitch * 2, coeff[-1], label = r'$\langle n_1\rangle^{{\infty}}={:.3g} \pm {:.3g}$'.format(coeff[-1], cov[-1]), s=fitmarkersize, marker='x', c='black')

        if searchkey('U', group[0]) != '5.0':
            ax2.plot( t_sec, f(t_sec, *coeff), label=r"$ \langle n_1(t)\rangle = A e^{-kt} + \langle n_1\rangle ^{\infty} $" + '\n' + r'$A = {:.3g}, k = {:.3g}$'.format(*coeff), linewidth=fitlinewidth, c='black')

    def fitcurrent2():

        def f(x, B, w, phi, A):
            return B * np.cos(w * (x + phi)) + A
        
        current = currents[0]
        time = times[0]
        sec1 = np.argwhere( time >= min(8, tswitch/2))
        sec2 = np.argwhere( time <= tswitch)

        idx = np.intersect1d(sec1, sec2)
        to_fit = current[idx]
        t_sec = time[idx]

        p0 = (0.01, 0.4, 6, 0.055)

        coeff, pcov = curve_fit(f, t_sec, to_fit, p0=p0)
        cov = np.sqrt(np.diag(pcov))

        print(cov)

        ax : plt.Axes =axes[ r + cat.index('current')][c]

        ax.scatter( tswitch , coeff[-1], label = r'$I_2^{{\infty}}={:.3g} \pm {:.3g}$'.format(coeff[-1], cov[-1]), s=fitmarkersize, marker='x', c='cyan')
        ax.plot( t_sec, f(t_sec, *coeff), label=r'$I_2(t) = \alpha  \cos [ \omega (t + \phi)] + I_2^{{\infty}}$' + '\n' + r'$\alpha = {:.3g}, \omega = {:.3g}, \phi = {:.3g}$'.format(*coeff[:-1]), linewidth=fitlinewidth, c= 'cyan')
        #ax.annotate( 'FIT STAGE 2 CURRENT : $  {:2f} * \cos ( {:2f} (t + {:2f})) + {:2f} $'.format(*coeff), (0.5, 0.2), color='black')

    def fitcurrent3():

        def f(x, B, w, phase, A):
            return B * np.cos(w * (x - CUTOFF + phase)) + A
        

        lens = [ len(c) if '34' in group[i] else 0 for i, c in enumerate(currents)]
        which = np.argsort(lens)[-1]

        p0 = (0.1, 0.5, 6, 0.06)

        if searchkey('U', group[0]) == '5.0':
            CUTOFF = tswitch * 1.5
        else:
            CUTOFF = tswitch * 1.1

        current = currents[which]
        time = times[which]
        idx = np.argwhere( time >= CUTOFF).flatten()

        to_fit = current[idx]
        t_sec = time[idx]
        

        coeff, pcov = curve_fit(f, t_sec, to_fit, p0=p0)
        cov = np.sqrt(np.diag(pcov))


        ax : plt.Axes =axes[ r + cat.index('current')][c]

        ax.scatter( tswitch , coeff[-1], label = r'$I_3^{{\infty}}={:.3g} \pm {:.3g}$'.format(coeff[-1], cov[-1]), s=fitmarkersize, marker='x', c='black')
        ax.plot( t_sec, f(t_sec, *coeff), label=r'$I_3(t) = \beta  \cos [ w (t + \delta)] + I_3^{{\infty}}$' + '\n' + r'$\beta = {:.3g}, w = {:.3g}, \delta = {:.3g}$'.format(*coeff[:-1]), c = 'black', linewidth=fitlinewidth)


    def plot_reference():
        Ls = set([searchkey('L', f) for f in group])
        Us = set([searchkey('U', f) for f in group])
        cnt = 0

        refgroup = [ 'L' + l + 'U' + u for l in Ls for u in Us]
        _, _, linestyles, markers= str_generator(refgroup, colorgroup, linestylegroup, markergroup, exclude_tag=exclude_tag)

        for L in Ls:
            for U in Us:
                
                linestyle = linestyles[cnt]
                marker = markers[cnt]

                if 'current' in cat:
                    plot_reference_current( L, U, linestyle, marker)

                if 'occ' in cat:
                    plot_reference_occ(L, U, linestyle, marker)

                cnt += 1

    def set_legend():

        shadecolor = 'grey'
        if 'current' in cat:
            cnt = cat.index('current')
            ax : plt.Axes =axes[ r + cnt][c]
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            
            #xprint(sorted(by_label.keys()))
            
            by_label = { key : by_label[key] for key in sorted(by_label.keys())}

            if fit:
                #ax.legend(by_label.values(), by_label.keys(), loc='lower center', bbox_to_anchor=(0.5, -0.2), fontsize=fontsize)
                ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)
            else:
                ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)

            if len(yscale_current):
                ax.set_ylim(yscale_current)
            ax.axvspan(tswitch  * 2- shadeoff, tswitch * 2, color=shadecolor, alpha=0.5)
            ax.tick_params(axis='both', labelsize=tickfontsize)

        if 'occ' in cat:

            if not current_occ_separate:
                cnt = cat.index('current')
                ax2 : plt.Axes = axes[ r + cnt ][c].twinx()

            else:
                cnt = cat.index('occ')
                ax2 : plt.Axes = axes[r+ cnt][c]

            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 1, box.height])

            # Put a legend to the right of the current axis

            handles, labels = ax2.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            by_label = { key : by_label[key] for key in sorted(by_label.keys())}
            if fit:
                #ax2.legend(by_label.values(), by_label.keys(), loc='lower center', bbox_to_anchor=(0.5, -0.2), fontsize=fontsize)
                ax2.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)
            else:
                ax2.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)

            ax2.axvspan(tswitch  * 2- shadeoff, tswitch * 2, color=shadecolor, alpha=0.5)
            ax2.tick_params(axis='both', labelsize=tickfontsize)

        
        if 'effE' in cat:
            cnt = cat.index('effE')
            ax : plt.Axes =axes[ r + cnt][c]

            ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)

            ax.axvspan(tswitch  * 2- shadeoff, tswitch * 2, color=shadecolor, alpha=0.5)
            ax.tick_params(axis='both', labelsize=tickfontsize)


    with PdfPages('plots/{}multi.pdf'.format(filetag)) as pdf:

        for page in inputfiles:

            cur_dict = inputfiles[page]
            yscale_current = cur_dict['yscale']
            cat :list = cur_dict['cat']

            if fit and 'effE' in cat:
                cat.remove('effE')
            groups : list[list[str]] = cur_dict['groups']
            group_ids = cur_dict['group_ids']
            exclude_tag : set = cur_dict['exclude_tag']

            # this defaults to none
            linestylegroup = cur_dict['linestylegroup']
            markergroup = cur_dict['markergroup']
            colorgroup = cur_dict['colorgroup']

            if len(groups) != len(group_ids):
                raise(ValueError('group length and ids dont match'))
            
            col = cur_dict['wrap']
            current_occ_separate = cur_dict['separate']
            multiplier = len(cat) - (current_occ_separate == False and ('occ' in cat and 'current' in cat) )

            row = ((len(groups) - 1) // col + 1) * multiplier

            if ('EE' in cat or 'bond' in cat) and max([ len(val) for val in groups]) > 1:
                raise(ValueError('more than 1 obj in each cat for EE'))
            
            if fit:
                width = 20
            else:
                width= 25 

            fig, axes = plt.subplots( row, col, figsize = (width * col, 4.5 * row))

            axes = [axes] if row == 1 else axes

            axes = [ [ax] for ax in axes] if  col == 1 else axes

            if current_occ_separate:
                currentcmap = mpl.cm.coolwarm

            else:
                currentcmap = mpl.cm.cool

            occcmap = mpl.cm.coolwarm
            markersize =15
            fitlinewidth = 10
            fitmarkersize = 300
            linewidth =5
            fontsize=35
            tickfontsize = 20

            # title location
            xtitle = 0.18
            ytitle = 0.5


            # we plot everything in the same group on one plot, plus multiplicity
            for i, group in enumerate(groups):
                
                

                group_id = group_ids[i]

                r =  (i // col) * multiplier
                c = i % col
                tswitch = float(searchkey('tswitch', group[0]))

                

                # we would change the context of all tags and labels, so we plot reference before
                ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(cat, group)
                tags, cvals, linestyles, markers= str_generator(group, colorgroup, linestylegroup, markergroup, exclude_tag=exclude_tag)
                

                if fit:

                    for func in (
                        fitocc,
                        fitcurrent2,
                        fitcurrent3
                    ):
                        
                        try:
                            func()
                        except RuntimeError:
                            print("what")

                    set_legend()

                if not fit:
                    plot_reference()
                


                for j, file in enumerate(group):
                    

                    time = times[j]
                    tag = tags[j]
                    cval = cvals[j]
                    linestyle = linestyles[j]
                    marker = markers[j]
                    
                    if 'bond' in cat:
                        plotbond()

                    if 'EE' in cat:
                        plotee()

                    if 'effE' in cat:
                        plotseff()

                    if 'current' in cat:
                        plotcurrent()
                        
                    if 'occ' in cat:
                        plotocc()


                if not fit:
                    set_legend()
                

                
            # for axr in axes:
            #     for a in axr:
            #         a.set_xticklabels([])
            #         a.set_yticklabels([])
            #         #a.set_aspect('equal')
                    
            #fig.suptitle(page)

            fig.tight_layout()
            fig.subplots_adjust(wspace=0, hspace=0)


            if save_individual:
                fig.savefig('plots/' + page + '.pdf')


            pdf.savefig(fig)




    return 0

def plotdensity(files):


    ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(['occ'], files)
    ts = [0.25]
    fig, axes = plt.subplots(len(files), len(ts), figsize=(30 * len(ts), 8 * len(files)))

    xtitle = 0.18
    ytitle = 0.5
    markersize = 300
    fontsize = 20
    for j, file in enumerate(files):
        occ = occs[j]

        tag = get_tag(file, exclude_tag={'ddposition', 'avg', 'L'})
        t = ts[0]
        ax2 : plt.Axes = axes[j]

        idx = np.argwhere( times[j] == t).flatten()[0]

        den = occ[idx, :-2]
        div = den.shape[0] // 2
        L = div - 2
        R = div + 2
        
        
        print(den)
        ax2.scatter(np.arange(1, L + 1 ), den[:L], label='Left', color='blue',
            marker='o', s=markersize
        )

        ax2.scatter(np.arange(L + 1, div + 1 ), den[L:div], label='LeftQPC', color='violet',
        marker='o', s=markersize
        )

        ax2.scatter(np.arange(div + 1, R + 1 ), den[div:R], label='RightQPC', color='orange',
            marker='o', s=markersize
        )

        ax2.scatter(np.arange(R + 1, div * 2 + 1 ), den[R:], label='Right', color='red',
        marker='o', s=markersize
        )

        ax2.set_ylabel(r'$\langle n_1\rangle$',  fontsize=fontsize)
        ax2.set_xlabel('Site number',  fontsize=fontsize)

        ax2.legend( fontsize=fontsize)
        ax2.tick_params(axis='both', labelsize=fontsize)
        ax2.set_title(r'$\langle n_k\rangle ({}),$'.format(t) + tag, fontsize=fontsize, x=xtitle, y=ytitle)

    fig.savefig('plots/compden.pdf')
# def comp_mf():


#     files = sorted(glob.glob("U*128*256*32"))

#     fig, ax = plt.subplots()

#     #refcur= np.loadtxt('meanfield/current')
#     #refocc = np.loadtxt('meanfield/occ')

#     ax : plt.Axes = ax
    
#     color = 'blue'
#     # ax.errorbar( refcur[:, 0], refcur[:, 1], yerr= [refcur[:, 2], refcur[:,3]], fmt='-o', color=color, label='Mean field current')

#     ax.set_xlabel('U')
#     ax.set_ylabel('Current (all)', color=color)
#     ax.tick_params(axis='y', labelcolor=color)

#     ax2 :plt.Axes = ax.twinx()

#     color = 'red'
#     # ax2.errorbar( refocc[:, 0], refocc[:, 1], yerr= [refocc[:, 2], refocc[:,3]], fmt='-x', color=color, label ='Mean field occ')

#     ax2.set_ylabel(r'$\langle n_1\rangle$', color=color)
#     ax2.set_title('Avg. current and lower dot occupation vs. U, Meanfield vs. MPS')
#     ax2.tick_params(axis='y', labelcolor=color)

#     cutoff = int(16/0.25)
#     rawcurrents = [   np.loadtxt(file + '/currentLR')[-cutoff:] for file in files]
#     currents = np.array([ [np.average(current), np.average(current) -np.amin(current), np.amax(current) - np.average(current)] for current in rawcurrents])
#     U = [float(file.split('_')[0][len('U'):]) for file in files]

#     rawocc = [   np.loadtxt(file + '/occ')[-cutoff:, -2] for file in files]
#     LT = np.array([ [np.average(occ), np.average(occ) -np.amin(occ), np.amax(occ) - np.average(occ)] for occ in rawocc])

        
#     ax.errorbar( U, currents[:, 0], yerr= [currents[:, 1], currents[:,2]], fmt='-o', color='purple', label='MPS current')
#     ax2.errorbar( U, LT[:, 0], yerr= [LT[:, 1], LT[:,2]], fmt='-x', color='orange', label='MPS occ')



#     ax.legend(loc ='lower left')
#     ax2.legend(loc='center right')

#     ax2.set_ylim(0, 1.1)
#     fig.savefig('plots/MF.pdf')




def Best_Fit32():

    save_individual = True
    extend_left = True
    shadeoff = 3
    yscale = np.array([-0.005, 0.15])

    exclude_tag = {'ddposition', 'avg', 'basis', 'dim'}

    ref = { 'linestylegroup' : { 'includeUFalse' : 'solid'},
                'colorgroup' : { 'mixedFalse': 0.2,
                                'includeUTrue' : 0.7,
                                'includeUFalse' : 0.95},
                'markergroup' : {'L64' : '',
                                'L34': ''}
    }

    search_dict = {

        'Fit' + 'Best32U2.0' : {'cat': ['occ', 'current'],
                                           'groups': [strsort('U2.0*L34*1024*includeUFalse*avgFalse*')], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Fit' +'Best32U3.0' : {'cat': ['occ', 'current'],
                                           'groups': [strsort('U3.0*L34*1024*includeUFalse*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Fit' + 'Best32U5.0' : {'cat': ['occ', 'current'],
                                           'groups': [strsort('U5.0*L32*1024*includeUFalse*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : [],
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           }
    

    }

    time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=True, save_individual=save_individual, extend_left=extend_left, shadeoff=shadeoff, MF=False)

def Best_Fit64():

    save_individual = True
    extend_left = True
    shadeoff = 5
    yscale = np.array([-0.005, 0.15])

    ref = { 'linestylegroup' : { 
                                'includeUFalse' : 'solid'},
                'colorgroup' : { 'mixedFalse': 0.2,
                                'includeUTrue' : 0.7,
                                'includeUFalse' : 0.95},
                'markergroup' : {'L64' : '',
                                'L34': ''}
    }
    exclude_tag  = {'ddposition', 'avg', 'basis', 'dim'}

    search_dict = {

        'Fit' + 'Best64U2.0' : {'cat': ['occ', 'current'],
                                           'groups': [strsort('U2.0*L64*512*includeUFalse*avgFalse*')], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Fit' +'Best64U3.0' : {'cat': ['occ', 'current'],
                                           'groups': [strsort('U3.0*L64*512*includeUFalse*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Fit' + 'Best64U5.0' : {'cat': ['occ', 'current'],
                                           'groups': [strsort('U5.0*L64*512*includeUFalse*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : [],
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           }
    

    }

    time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=True, save_individual=save_individual, extend_left=extend_left, shadeoff=shadeoff, MF=False)



def Best_32wrap(MF=False):


    save_individual = True
    extend_left = True
    shadeoff = 3
    yscale = np.array([-0.005, 0.15])

    cat = ['occ', 'current'] + (['effE'] if not MF else [])

    ref = { 'linestylegroup' : { 'mixedFalse': 'solid',
                                'includeUTrue' : 'dotted',
                                'includeUFalse' : 'dashed'},
                'colorgroup' : { 'mixedFalse': 0.2,
                                'includeUTrue' : 0.7,
                                'includeUFalse' : 0.95},
                'markergroup' : {'L32' : ''}
    }

    exclude_tag = {'ddposition', 'avg', 'L', 'dim'}

    search_dict = {

        'Best32U2.0' : {'cat': cat,
                                           'groups': [strsort('U2.0*L32*1024*avgFalse*')], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Best32U3.0' : {'cat': cat,
                                           'groups': [strsort('U3.0*L32*1024*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Best32U5.0' : {'cat': cat,
                                           'groups': [strsort('U5.0*L32*1024*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : [],
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           }
    

    }

    time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=False, save_individual=save_individual, extend_left=extend_left, shadeoff=shadeoff, MF=MF)


def Best_64wrap(MF=False):

    cat = ['occ', 'current'] + (['effE'] if not MF else [])
    save_individual = True
    extend_left = True
    shadeoff = 5
    yscale = np.array([-0.005, 0.15])

    ref = { 'linestylegroup' : { 'mixedFalse': 'solid',
                                'includeUTrue' : 'dotted',
                                'includeUFalse' : 'dashed'},
                'colorgroup' : { 'mixedFalse': 0.2,
                                'includeUTrue' : 0.7,
                                'includeUFalse' : 0.95},
                'markergroup' : {'L32' : ''}
    }

    exclude_tag = {'ddposition', 'avg', 'L', 'dim'}

    search_dict = {

        'Best64U.0' : {'cat': cat,
                                           'groups': [strsort('U2.0*L64*512*avgFalse*')], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                            'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Best64U3.0' : {'cat': cat,
                                           'groups': [strsort('U3.0*L64*512*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'Best64U5.0' : {'cat': cat,
                                           'groups': [strsort('U5.0*L64*512*avgFalse*') ], 
                                           'group_ids': [''], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : [],
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           }
    

    }

    time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=False, save_individual=save_individual, extend_left=extend_left, shadeoff=shadeoff, MF=MF)


def Best_3234wrap(MF=True):

    filetag = 'BestwMF' + str(MF)
    fit = False
    save_individual = True
    extend_left = True
    shadeoff = 3

    cat = ['occ', 'current'] + (['effE'] if not MF else [])

    yscale = np.array([-0.005, 0.15])

    ref = { 'linestylegroup' : { 'L32': 'solid',
                                'L34' : 'dashed',
                                },
                'colorgroup' : { 'L32': 0.1,
                                'L34' : 0.9},
                'markergroup' : {'L32' : '',
                                'L34': ''}
    }

    exclude_tag = {'ddposition', 'avg', 'basis'}

    search_dict = {


    'BestwMF' +  str(MF) + '2.0' : {'cat': cat,
                                        'groups': [strsort('U2.0*L3*1024*includeUFalse*avgFalse*')], 
                                        'group_ids': [''], 
                                        'separate' : True,
                                        'wrap' : 1,
                                        'yscale' : yscale,
                                        'exclude_tag' : exclude_tag,
                                        'linestylegroup' : ref['linestylegroup'],
                                        'colorgroup' : ref['colorgroup'],
                                        'markergroup' : ref['markergroup']
                                        },
    'BestwMF' + str(MF) + '3.0' : {'cat': cat,
                                        'groups': [strsort('U3.0*L3*1024*includeUFalse*avgFalse*') ], 
                                        'group_ids': [''], 
                                        'separate' : True,
                                        'wrap' : 1,
                                        'yscale' : yscale,
                                        'exclude_tag' : exclude_tag,
                                        'linestylegroup' : ref['linestylegroup'],
                                        'colorgroup' : ref['colorgroup'],
                                        'markergroup' : ref['markergroup']
                                        },
    'BestwMF' + str(MF)+ '5.0' : {'cat': cat,
                                        'groups': [strsort('U5.0*L3*1024*includeUFalse*avgFalse*') ], 
                                        'group_ids': [''], 
                                        'separate' : True,
                                        'wrap' : 1,
                                        'yscale' : [],
                                        'exclude_tag' : exclude_tag,
                                        'linestylegroup' : ref['linestylegroup'],
                                        'colorgroup' : ref['colorgroup'],
                                        'markergroup' : ref['markergroup']
                                        }
    }

    time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=fit, save_individual=save_individual, extend_left=extend_left, shadeoff=shadeoff, MF=MF)



def test_0125(MF=True):

    filetag = 'test0125withEFF' 
    save_individual = True
    extend_left = True
    shadeoff = 3
    yscale = np.array([-0.005, 0.15])

    cat = ['occ', 'current', 'effE'] #+ (['effE'] if not MF else [])

    ref = { 'linestylegroup' : { 'mixedFalse': 'solid',
                                'includeUTrue' : 'dotted',
                                'includeUFalse' : 'dashed'},
                'colorgroup' : { 'mixedFalse': 0.2,
                                'includeUTrue' : 0.7,
                                'includeUFalse' : 0.95},
                'markergroup' : {'1024' : 'o',
                                 '512' : ''}
    }

    exclude_tag = {'ddposition', 'avg'}

    search_dict = {

        'BestU0.5' : {'cat': cat,
                                           #'groups': [strsort('U0.5*34*512*')], 
                                           'groups': [strsort('U0.5*34*')], 
                                           'group_ids': ['$U=0.5$'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'BestU1.0' : {'cat': cat,
                                           #'groups': [strsort('U1.0*34*512*') ], 
                                           'groups': [strsort('U1.0*34*')], 
                                           'group_ids': ['$U=1.0$'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'BestU1.5' : {'cat': cat,
                                           #'groups': [strsort('U1.5*34*512*') ], 
                                           'groups': [strsort('U1.5*34*') ], 
                                           'group_ids': ['$U=1.5$'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : [],
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },

        'BestU0.5L66' : {'cat': cat,
                                           #'groups': [strsort('U0.5*34*512*')], 
                                           'groups': [strsort('U0.5*66*')], 
                                           'group_ids': ['$U=0.5$'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'BestU1.0L66' : {'cat': cat,
                                           #'groups': [strsort('U1.0*34*512*') ], 
                                           'groups': [strsort('U1.0*66*')], 
                                           'group_ids': ['$U=1.0$'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : yscale,
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           },
        'BestU1.5L66' : {'cat': cat,
                                           #'groups': [strsort('U1.5*34*512*') ], 
                                           'groups': [strsort('U1.5*66*') ], 
                                           'group_ids': ['$U=1.5$'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'yscale' : [],
                                           'exclude_tag' : exclude_tag,
                                           'linestylegroup' : ref['linestylegroup'],
                                           'colorgroup' : ref['colorgroup'],
                                            'markergroup' : ref['markergroup']
                                           }
    

    }

    time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=False, save_individual=save_individual, extend_left=extend_left, shadeoff=shadeoff, MF=MF)


def test_den():


    files = strsort('U2.0*L64*512*mixedTrue*avgFalse*')
    plotdensity(files=files)



if __name__ == '__main__':


    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 

    #filetag = 'FIT64'
    #filetag = 'L64dim512'
    filetag = 'BestwMF'
    fit = False
    save_individual = True
    extend_left = True
    shade = True
    MF = True

    #sample : U2.0_L32_TEdim1024_tswitch16.0_includeUTrue_orderingSORTED_mixedTrue_avgFalse_ddpositionR

    #test_den()

    # Best_32wrap()
    # Best_3234wrap(MF=False)
    # Best_Fit32()
    # Best_64wrap()
    # Best_Fit64()
    # Best_3234wrap(MF=True)

    test_0125()

    #current_diff(inputoverride=strsort('*2.0*L64*'))
    #time_bond_ee(inputfiles=search_dict, filetag=filetag, fit=fit, save_individual=save_individual, extend_left=extend_left, shade=shade, MF=MF)


