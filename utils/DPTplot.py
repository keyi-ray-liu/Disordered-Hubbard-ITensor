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


def str_generator(files, linestylegroup, markergroup):

    tags = sorted([ get_tag(file) for file in files])
    cval = np.linspace(0.1, 0.9, len(tags))

    linestyle_str = [
     'solid',      # Same as (0, ()) or '-'
     'solid',    # Same as (0, (1, 1)) or ':'
     'dashed',    # Same as '--'
     'dashdot']  # Same as '-.'
    
    marker_str = ['o', 'x', '^']

    comb = product(linestyle_str, marker_str)


    if linestylegroup == {}:
        linestyles = [ val[0] for val in comb]

    else:

        linestyles = []

        for file in files:
            for key in linestylegroup:
                if key in file:
                    linestyles += [ linestylegroup[key]]


    if markergroup == {}:
        markers = [ val[1] for val in comb] 

    else:

        markers = []

        for file in files:
            for key in markergroup:
                if key in file:
                    markers += [ markergroup[key]]


    return tags, cval, linestyles, markers


def get_tag(file):

    tag = {}

    if 'mix' not in file or 'mixedFalse' in file:
        tag['basis'] = 'spatial'
        tag['includeQPC'] = 'T'
        

    else: 
        if 'random' in file:
            tag['basis'] = 'mixed random'

        elif 'orderingABSSORTED' in file:
            tag['basis'] = 'mixed abs ordered'

        #elif 'orderingSORTED' in file:
        else:
            tag['basis'] = 'mixed ordered'

        if 'includeUFalse' in file:
            tag['includeQPC'] = 'F'

        else:
            tag['includeQPC'] = 'T'

    if 'ddpositionL' in file:
        tag['DDpos'] = 'L'

    elif 'ddpositionR' in file:
        tag['DDpos'] = 'R'

    elif 'ddpositionM' in file:
        tag['DDpos'] = 'M'

    if 'avgTrue' in file:
        tag['avg'] = 'T'

    else:
        tag['avg'] = 'F'



    tag['L'] = searchkey('L', file)
    res = ', '.join([ key + '=' + tag[key] for key in tag])
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
        tag, linestyle = get_tag(file)

        U = searchkey('U', file)

        tag = tag +  'bd = ' + ''.join([ searchkey('TEdim', f) for f in file.split('_')])
        
        
        curdim = searchkey('TEdim', file)

        
        replaced = file.replace(curdim, maxdim)
        ref = np.loadtxt( replaced + '/currentLR')

        if 'abs' not in tag and '32' not in file:
            ref = np.concatenate( ([0.0 for _ in range(12)], ref))

        time = times[i][:ref.shape[0]]

        cval = set_color(files, curdim, 'TEdim')
        

        if replaced != file:
            current = currents[i]

            ax : plt.Axes = axes[axdict[U]][0]

            div = np.where( ref > 1e-6, ref, 1)
            diff = (current[:ref.shape[0]] - ref)/div
            ax.plot( time, diff, label=tag, c=cmap(cval), linestyle=linestyle)

            ax: plt.Axes = axes[axdict[U]][1]

            diff = np.abs(diff)
            ax.plot( time, diff, label=tag, c=cmap(cval), linestyle=linestyle)



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
    return fig




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




def time_bond_ee( inputfiles : dict = {},  filetag = 'multi' ):


    def plotcurrent():

        cnt = cat.index('current')
        current = currents[j]
        ax : plt.Axes =axes[ r + cnt][c]
        color = currentcmap(cval)

        if time[0] > -tswitch:
            ctime = np.concatenate(( [-tswitch], time))
            current = np.concatenate( ([0], current))
        
        if not current_occ_separate:
            axiscolor = color
            ax.tick_params(axis='y', labelcolor=color)

        else:
            axiscolor='black'

        ax.plot(ctime, current, label=r'$I_{L\rightarrow R}$ , MPS ' + tag, color=color, 
            linestyle=linestyle,
            marker=marker, markersize=markersize
            )
        
        ax.vlines([0, tswitch], np.zeros(2), 1.1 * np.amax(current) * np.ones(2), linestyles='dotted')
        ax.set_xlabel('time $(1/\omega_0)$',  fontsize=fontsize)
        ax.set_ylabel('Current (all)', color=axiscolor,  fontsize=fontsize)

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        by_label = { key : by_label[key] for key in sorted(by_label.keys())}
        ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)

        ax.set_title('Current vs. t. {}'.format(group_id),  fontsize=fontsize)

    def plot_reference_current(L, U):

        
        ax : plt.Axes =axes[ r + cat.index('current')][c]
        refcolor='green'
        ref = np.loadtxt('/Users/knl20/Desktop/Code/TN/non-interacting/{}currentCC{}'.format(L, U))

        reftime = np.arange(0, tswitch, 0.25)
        ax.plot( reftime, ref[:reftime.shape[0]], label=r'$I_{L\rightarrow R}$, Stage 2, non-interacting exact, '+ '$L={}$'.format(L), color=refcolor,
                linestyle=linestyle, markersize=markersize, marker=marker
                )

        refcolor='violet'
        refraw = np.loadtxt('/Users/knl20/Desktop/Code/TN/LT/meanfield/Current_dpt_mf_N{}_mu0.5_U{}_tswitch{}'.format(L, U, tswitch))

        reftime = refraw[:, 0]
        ref = refraw[:, 3]

        ax.plot( reftime, ref[:reftime.shape[0]], label=r'$I_{L\rightarrow R}$, MF, ' + '$L={}$'.format(L), color=refcolor,
                linestyle=linestyle, marker=marker, markersize=markersize
                )

    def plot_reference_occ(L, U):

        if not current_occ_separate:
            ax2 : plt.Axes = axes[ r + cat.index('current')][c].twinx()

        else:
            ax2 : plt.Axes = axes[r+ cat.index('occ')][c]

        refcolor='violet'
        refraw = np.loadtxt('/Users/knl20/Desktop/Code/TN/LT/meanfield/Current_dpt_mf_N{}_mu0.5_U{}_tswitch{}'.format(L, U, tswitch))

        reftime = refraw[:, 0]
        ref = refraw[:, 5]

        ax2.plot( reftime, ref[:reftime.shape[0]], label='lower dot density, MF, L={}'.format(L), color=refcolor,
                linestyle=linestyle, marker=marker, markersize=markersize
                )

    def plotocc():
        occ = occs[j]
        color = occcmap(cval)

        if not current_occ_separate:
            ax2 : plt.Axes = axes[ r + cat.index('current')][c].twinx()

        else:
            ax2 : plt.Axes = axes[r+ cat.index('occ')][c]
        

        if 'ddpositionL' in file:
            ddlow = 0

        elif 'ddpositionM' in file or 'avgTrue' in file: 
            ddlow = occ.shape[-1]//2 - 1

        else:
            ddlow = -2

        LT1 = occ[:, ddlow]
        obj,  = ax2.plot(time, LT1, label='lower dot density, ' + tag, color=color,
        linestyle=linestyle, marker=marker, markersize=markersize
        )

        if not current_occ_separate:
            axiscolor = color
            ax2.tick_params(axis='y', labelcolor=color)

        else:
            axiscolor='black'

        
        ax2.set_ylabel(r'$\langle n_S\rangle$', color=axiscolor,  fontsize=fontsize)
        ax2.set_xlabel('time $(1/\omega_0)$',  fontsize=fontsize)
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 1, box.height])

        # Put a legend to the right of the current axis

        handles, labels = ax2.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        by_label = { key : by_label[key] for key in sorted(by_label.keys())}
        ax2.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)

        ax2.vlines([0, tswitch], np.zeros(2), 1.1 * np.ones(2), linestyles='dotted')
        ax2.set_title('Lower Dot Occupation vs. t. {}'.format(group_id), fontsize=fontsize)

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

        if time[0] > -tswitch:
            ctime = np.concatenate(( [-tswitch], time))
            seff = np.concatenate( ([seff[0]], seff))


        ax.plot(ctime, seff, label='Seff, MPS' + tag, color=color, 
            linestyle=linestyle,
            marker=marker, markersize=markersize
            #marker='x', s=5
            )
        
        ax.vlines([0, tswitch], np.zeros(2), 1.1 * np.amax(seff) * np.ones(2), linestyles='dotted')
        ax.set_xlabel('time $(1/\omega_0)$',  fontsize=fontsize)
        ax.set_ylabel('$S_{eff}$',  fontsize=fontsize)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.6), fontsize=fontsize)
        ax.set_title('Effectly Entropy vs. t, {}'.format(group_id), fontsize=fontsize)

    # input files are grouped by dicts, 
    # each dict: 
    # input = { "page1 title": {"group 1 on page1" : [], "group 2 on page1 : []", etc.}
    #


    with PdfPages('plots/{}multi.pdf'.format(filetag)) as pdf:

        for page in inputfiles:

            cur_dict = inputfiles[page]

            cat :list = cur_dict['cat']
            groups : list[list[str]] = cur_dict['groups']
            group_ids = cur_dict['group_ids']

            # this defaults to none
            linestylegroup = cur_dict['linestylegroup']
            markergroup = cur_dict['markergroup']

            if len(groups) != len(group_ids):
                raise(ValueError('group length and ids dont match'))
            
            col = cur_dict['wrap']
            current_occ_separate = cur_dict['separate']
            multiplier = len(cat) - (current_occ_separate == False and ('occ' in cat and 'current' in cat) )

            row = ((len(groups) - 1) // col + 1) * multiplier
            
            if 'EE' in cat or 'bond' in cat and max([ len(val) for val in groups]) > 1:
                raise(ValueError('more than 1 obj in each cat for EE'))
            
            fig, axes = plt.subplots( row, col, figsize = (30 * col, 7 * row))

            axes = [ [ax] for ax in axes] if  col == 1 else axes

            if current_occ_separate:
                currentcmap = mpl.cm.coolwarm

            else:
                currentcmap = mpl.cm.cool

            occcmap = mpl.cm.coolwarm
            markersize = 5
            fontsize=20

            # we plot everything in the same group on one plot, plus multiplicity
            for i, group in enumerate(groups):

                group_id = group_ids[i]

                r =  (i // col) * multiplier
                c = i % col

                ees, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax = dataloader(cat, group)
                tags, cvals, linestyles, markers= str_generator(group, linestylegroup, markergroup)

                tswitch = float(searchkey('tswitch', group[0]))
                

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
                        plot_reference_current( searchkey('L', file), searchkey('U', file))
                        plotcurrent()
                        
                    if 'occ' in cat:
                        plot_reference_occ( searchkey('L', file), searchkey('U', file))
                        plotocc()

                

            fig.suptitle(page)
            fig.tight_layout()
            pdf.savefig(fig)




    return 0



def comp_mf():


    files = sorted(glob.glob("U*128*256*32"))

    fig, ax = plt.subplots()

    #refcur= np.loadtxt('meanfield/current')
    #refocc = np.loadtxt('meanfield/occ')

    ax : plt.Axes = ax
    
    color = 'blue'
    # ax.errorbar( refcur[:, 0], refcur[:, 1], yerr= [refcur[:, 2], refcur[:,3]], fmt='-o', color=color, label='Mean field current')

    ax.set_xlabel('U')
    ax.set_ylabel('Current (all)', color=color)
    ax.tick_params(axis='y', labelcolor=color)

    ax2 :plt.Axes = ax.twinx()

    color = 'red'
    # ax2.errorbar( refocc[:, 0], refocc[:, 1], yerr= [refocc[:, 2], refocc[:,3]], fmt='-x', color=color, label ='Mean field occ')

    ax2.set_ylabel(r'$\langle n_S\rangle$', color=color)
    ax2.set_title('Avg. current and lower dot occupation vs. U, Meanfield vs. MPS')
    ax2.tick_params(axis='y', labelcolor=color)

    cutoff = int(16/0.25)
    rawcurrents = [   np.loadtxt(file + '/currentLR')[-cutoff:] for file in files]
    currents = np.array([ [np.average(current), np.average(current) -np.amin(current), np.amax(current) - np.average(current)] for current in rawcurrents])
    U = [float(file.split('_')[0][len('U'):]) for file in files]

    rawocc = [   np.loadtxt(file + '/occ')[-cutoff:, -2] for file in files]
    LT = np.array([ [np.average(occ), np.average(occ) -np.amin(occ), np.amax(occ) - np.average(occ)] for occ in rawocc])

        
    ax.errorbar( U, currents[:, 0], yerr= [currents[:, 1], currents[:,2]], fmt='-o', color='purple', label='MPS current')
    ax2.errorbar( U, LT[:, 0], yerr= [LT[:, 1], LT[:,2]], fmt='-x', color='orange', label='MPS occ')



    ax.legend(loc ='lower left')
    ax2.legend(loc='center right')

    ax2.set_ylim(0, 1.1)
    fig.savefig('plots/MF.pdf')


if __name__ == '__main__':


    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 


    search_dict = {
        'Compare L34 and L32, dim=1024, U=2.0' : {'cat': ['occ', 'current', 'effE'],
                                           'groups': [strsort('*2.0*1024*')], 
                                           'group_ids': ['U = 2.0'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'linestylegroup' : { '34': 'solid',
                                                               '32' : 'dashed'},
                                            'markergroup' : {'34': 'o',
                                                             '32': 'x'}
                                           },
        'Compare L34 and L32, dim=1024, U=3.0' : {'cat': ['occ', 'current', 'effE'],
                                           'groups': [strsort('*3.0*1024*') ], 
                                           'group_ids': ['U = 3.0'], 
                                           'separate' : True,
                                           'wrap' : 1,
                                           'linestylegroup' : { '34': 'solid',
                                                               '32' : 'dashed'},
                                            'markergroup' : {'34': 'o',
                                                             '32': 'x'}
                                           },

    }




    time_bond_ee(inputfiles=search_dict, filetag='comp1024')


