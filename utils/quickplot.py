from ntpath import realpath
from sys import platform
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib import cm
from scipy.optimize import curve_fit
import sys
import matplotlib.backends.backend_pdf
import os
import glob
from matplotlib.ticker import MaxNLocator

def syspre():
    if platform == 'win32':
        syspre = 'C:/Users/Ray/iCloudDrive'

    else:
        syspre = '/Users/rayliu'

    return syspre

def lengthcompplot():

    ed = defaultdict()
    mps = defaultdict()

    ed[12] = -19.309357049818
    mps[12] = -19.309357048327

    ed[14] = -23.518369349887
    mps[14] = -23.518369346716

    ed[16] = -27.866218760560
    mps[16] = -27.866218756627

    ed[18] = -32.335906388135
    mps[18] = -32.335906382640

    ed[20] = -36.914153564647
    mps[20] = -36.914153558694

    ed[22] =  -41.590298379263
    mps[22] = -41.590298372306

    ed[24] = -46.35559011
    mps[24] = -46.355590102471

    ed[26] = -51.202716787557
    mps[26] = -51.202716777823

    fig, ax = plt.subplots()

    ax.set_xlabel('Length of the chain')
    ax.set_title('Comparison of MPS vs. ED GS energies at half-filling')
    color= 'tab:red'
    key = sorted(ed.keys())
    abserr = [abs( ed[k] - mps[k]) for k in key]
    relerr = [abs( ed[k] - mps[k])/ abs(ed[k]) for k in key]
    ax.set_ylabel('Absolute Error', color=color)
    ax.scatter( key, abserr, label='Abs. err', color=color)
    ax.set_ylim( min(abserr) * 0.9, max(abserr) * 1.1)
    ax.tick_params(axis= 'y', labelcolor=color)

    ax2 = ax.twinx()
    color ='tab:blue'
    ax2.set_ylabel('Relative Error', color=color)
    ax2.scatter( key, relerr, label= 'Rel. err', color=color)
    ax2.set_ylim( min(relerr) * 0.9, max(relerr) * 1.1)
    ax2.tick_params(axis= 'y', labelcolor=color)

    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.show()


    sweep = np.arange(2, 102, 2)
    energies = [-40.889374506862275, -41.513276484797316, -41.58386147065303, -41.58919292318066, -41.58980726872991, -41.59024313968609, -41.59027647602754, -41.590289904686244, -41.59029483388616, -41.590297078841424, -41.59029768381566, -41.590298124041375, -41.59029825481698, -41.59029830719316, -41.590298352746494, -41.590298362538874, -41.59029837032205, -41.5902983715727, -41.59029837230328, -41.59029837230632, -41.59029837230635, -41.590298372306385, -41.590298372306215, -41.59029837230622, -41.59029837230629, -41.59029837230616, -41.59029837230621, -41.59029837230626, -41.59029837230625, -41.590298372306314, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628, -41.59029837230628]

    fig, ax = plt.subplots()
    ax.set_title('Difference in energy vs max bond dimension, L=22')
    ax.set_ylabel('Difference between ED and DMRG energy')
    ax.set_xlabel('Max bond dim')
    #ax.scatter(sweep, energies, s=5, label='MPS')
    #ax.scatter(sweep, [-41.590298379263] * len(sweep), c='r', s=5, label='ED')

    ax.scatter(sweep, abs( np.array(energies) - np.array([-41.590298379263] * len(sweep))))
    plt.yscale('log')
    plt.show()

    syssize = np.arange(12, 102, 2)
    maxsweep = [24.0, 30.0, 37.0, 43.0, 46.0, 60.0, 57.0, 71.0, 65.0, 86.0, 80.0, 78.0, 81.0, 75.0, 71.0, 71.0, 70.0, 78.0, 82.0, 90.0, 76.0, 74.0, 70.0, 79.0, 77.0, 114.0, 80.0, 84.0, 91.0, 87.0, 211.0, 88.0, 92.0, 90.0, 107.0, 171.0, 89.0, 114.0, 202.0, 118.0, 95.0, 176.0, 420.0, 217.0, 257.0]


    fig, ax = plt.subplots()
    ax.set_title('Max bond dim for energy convergence vs. system size, 1D, extended Hubbard, half-filling')
    ax.set_ylabel('Max bond dim')
    ax.set_xlabel('System size')
    ax.scatter( syssize, maxsweep)
    plt.show()

def weightplot():
    errs = {}
    #tag = ['1', '10', '100', 'dyna']
    tag = ['1', '10', '100']

    #m = ['x', 'v', '^', '+']
    m = ['x', 'v', '^']
    individual = 0
    

    for lam in range(len(tag)):
        fig, ax = plt.subplots()
        if platform == 'darwin':
            prefix = '/Users/rayliu/Desktop/Code/DisorderedML/RayGenerator/testnewgen/zero/'
            mpsprefix = '/Users/rayliu/Desktop/Code/iTensor/collect/weighttest/26weight{}/'.format( tag[lam])

        
        #ed = np.loadtxt( prefix + 'energy')[0]
        ed = np.loadtxt( prefix + '26energy')

        exMPS = np.loadtxt( mpsprefix + 'ex')

        
        ax.set_xlabel('Excited states')
        ax.set_title('Comparison of MPS vs. ED Excited states search, L=26, $\lambda$ = {}'.format(tag[lam]))
        color= 'tab:red'
        abserr = [abs( exMPS[i] - ed[i]) for i in range(len(exMPS))]
        relerr = [abs( exMPS[i] - ed[i])/ abs(ed[i]) for i in range(len(exMPS))]

        errs[lam] = relerr
        ref = np.arange(1, len(exMPS) + 1, 1)
        ax.set_ylabel('Absolute Error', color=color)
        ax.scatter( ref, abserr, label='Abs err', color=color)
        #ax.set_ylim( min(abserr) * 0.9, max(abserr) * 1.1)
        ax.tick_params(axis= 'y', labelcolor=color)
        ax.set_yscale('log')

        ax2 = ax.twinx()
        color ='tab:blue'
        ax2.set_ylabel('Relative Error', color=color)
        ax2.scatter( ref, relerr, label= 'Rel err', color=color)
        ax2.set_ylim( min(relerr) * 0.9, max(relerr) * 1.1)
        ax2.tick_params(axis= 'y', labelcolor=color)
        ax2.set_yscale('log')

        if individual:
            plt.show()



    plt.close("all")
    fig, ax = plt.subplots()
    for lam in errs:

        ax.scatter( np.arange(1, len(errs[lam]) + 1), errs[lam], label='$\lambda$= {}'.format(tag[lam]), marker=m[lam])

    #sub = np.loadtxt('/Users/rayliu/Desktop/Code/MPS/obs/ranges/26full/exsubtract')
    #suberr = [abs( sub[i] - ed[i])/ abs(ed[i]) for i in range(len(sub))]
    #ax.scatter( np.arange(1, len(sub) +1), suberr, label='subtraction based, openmps' )
    ax.set_title('Comparison of relative errors, L=26')
    ax.set_ylabel('Relative Error')
    ax.set_xlabel('Excited state')
    ax.set_yscale('log')
    plt.legend()
    plt.show()


def plot2d():

    if platform == 'darwin':
        edprefix = '/Users/rayliu/Desktop/Code/PEPS/etED/comp/'
        mpsprefix = '/Users/rayliu/Desktop/Code/iTensor/collect/'

    fig, ax= plt.subplots()
    ed = np.loadtxt( edprefix + 'energy' )
    mps = np.loadtxt( mpsprefix  + '3x3weight10/ex')

    relerr = [abs( mps[i] - ed[i])/ abs(ed[i]) for i in range(len(mps))]

    ax.scatter ( np.arange(1, len(mps) + 1), relerr)
    ax.set_title( 'Comparison of difference between ED and MPS, 3x3, N=4')
    ax.set_ylabel('Relative Error')
    ax.set_xlabel('Eigenstate')
    ax.set_ylim( max(min(relerr), 1e-10) , max(relerr))

    ax.set_yscale('log')
    
    plt.show()

def longplot():
    

    if platform == 'darwin':
        prefix = '/Users/rayliu/Desktop/Code/iTensor/collect/'

    else:
        prefix = 'C:/Users/Ray/iCloudDrive/Desktop/Code/iTensor/collect/'

    cmap = cm.get_cmap('coolwarm')

    for L in (40, 60):

        fig, ax = plt.subplots()

        energy = np.loadtxt( prefix + '{}weight10/ex'.format(L))
        ref = np.arange(1, len(energy) + 1)
        #ax.scatter(  ref, energy, label = 'Energy')

        for i in range(len(energy)):
 
            val = i / (len(energy) - 1)
            ax.scatter( energy[i], energy[i] - energy[0], label = 'Energy difference to GS', color=cmap(val))
        
        ax.set_title('Excited state energy level, L = {} (Color indicates the order a state appears in the search)'.format(L))
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Energy difference to GS (eV)')

    plt.show()

def plasmonplot():

    def line(x, a, b):
        return a * x + b

    if platform == 'darwin':
        prefix = '/Users/rayliu/Desktop/Code/iTensor/collect/plasmon/ee1.0/'

    else:
        prefix = 'C:/Users/Ray/iCloudDrive/Desktop/Code/iTensor/collect/plasmon/ee1.0/'

    begin = 10
    end = 80


    plasmon = np.loadtxt(prefix + 'plasmonenergy')

    ref = np.arange(begin, end + 1, 2)

    logref, logplsm = np.log10(ref), np.log10(plasmon)
    fit, cov = curve_fit(line, logref[-5:], logplsm[-5:])
    print(fit)
    fig, ax = plt.subplots(3, figsize=(5, 8))

    ax[0].set_title('Plasmon plots, $\lambda_{ee} = \lambda_{ne} = 1.0$')

    ax[0].scatter(ref, plasmon)
    ax[0].set_xlabel('L')
    ax[0].set_ylabel('Plasmon energy')

    ax[1].scatter(1/ref, plasmon)
    ax[1].set_xlabel('1/L')
    ax[1].set_ylabel('Plasmon energy')

    ax[2].scatter( logref, logplsm)

    a, b= fit
    lineref = np.linspace( min(logref), max(logref), 100)
    ax[2].plot( lineref, line(lineref, a, b), c='r' , label = ' y = {:.3f}x + {:.3f}'.format(a, b))
    ax[2].legend()


    ax[2].set_xlabel('$log_{10} \  L$')
    ax[2].set_ylabel('$log_{10} \  E_{plasmon}$')

    fig.tight_layout()
    plt.show()

def highplasmonplot():

    lam = sys.argv[1]
    header = np.loadtxt('/Users/rayliu/Desktop/Code/iTensor/collect/plasmon/high' + lam + '/header')
    raw = np.loadtxt('/Users/rayliu/Desktop/Code/iTensor/collect/plasmon/high' + lam + '/energy')

    new = np.zeros( raw.shape)
    for i, l in enumerate(raw):
        new[i] = (l - l[0]) / ( l[1] - l[0])


    else:

        if platform == 'win32':
            syspre = 'C:/Users/Ray/iCloudDrive'

        else:
            syspre = '/Users/rayliu'

        plotloc = syspre + '/Desktop/Code/iTensor/plots/'

        fig, ax = plt.subplots()
        ax.scatter( header, new[:, 3], facecolors='none', edgecolors='r')
        ax.scatter( header, new[:, 6], facecolors='none', edgecolors='b')
        ax.scatter( header, new[:, 11], facecolors='none', edgecolors='orange' )
        ax.scatter( header, new[:, 19], facecolors='none', edgecolors='black' )

        ref = np.linspace( 0, max(header) + 5, 100)
        ax.plot(ref, np.ones( len(ref)) * 2, '--', c='r')
        ax.plot(ref, np.ones( len(ref)) * 3, '--', c='b')
        ax.plot(ref, np.ones( len(ref)) * 4, '--', c='orange')
        ax.plot(ref, np.ones( len(ref)) * 5, '--', c='black')

        ax.set_xlabel('L')
        ax.set_ylabel('Plasmon number')
        ax.set_title('Plasmon number vs. length of chain, $\lambda =$ {}'.format(lam))
        #plt.show()
        plt.savefig(plotloc + 'ex_plasmon_lam' + lam + '.png', dpi=600)

        fig, ax = plt.subplots()
        ax.scatter( header, new[:, 3] - 2, c='r', marker='o', label= '$N_{pl} = 2$' )
        ax.scatter( header, new[:, 6] - 3, c='b', marker = '^', label= '$N_{pl} = 3$' )
        ax.scatter( header, new[:, 11] - 4, c='orange', marker = 'x' , label= '$N_{pl} = 4$' )
        ax.scatter( header, new[:, 19] - 5, c='black' , marker = 'v', label= '$N_{pl} = 5$' )
        ax.legend()



        ax.set_xlabel('L')
        ax.set_ylabel('$\Delta N_{plasmon}$')
        ax.set_title('Difference in plasmon number to integer, vs. length of chain, $\lambda =$ {}'.format(lam))
        plt.savefig(plotloc + 'ex_plasmon__diff_lam' + lam +'.png' , dpi=600)
        #plt.show()

def extrendplot():

    pdf = matplotlib.backends.backend_pdf.PdfPages( syspre() + '/Desktop/Code/iTensor/plots/trendplot.pdf' )
    plasmondir = syspre () + '/Desktop/Code/iTensor/collect/plasmon/'

    
    for dir in sorted( glob.glob( plasmondir  + 'high*'), key= lambda x: float(x[ len(plasmondir) + 4:])):

        lam = dir[ len(plasmondir) + 4:]
        header = np.loadtxt(  dir + '/header')
        raw = np.loadtxt( dir + '/energy')

        new = np.zeros( raw.shape)
        for i, l in enumerate(raw):
            new[i] = (l - l[0]) / ( l[1] - l[0])

        fig, ax = plt.subplots()
        for ex in range(new.shape[1]):
            ax.plot( header, new[:, ex])

        ref = np.linspace( 5, 100, 5)
        for num in range(1, int(np.amax(new)) + 1):
            ax.plot(ref, num * np.ones(len(ref)), ':', c='black')

        ax.set_xlabel('Length of chain')
        ax.set_ylabel('Plasmon number')
        ax.set_xlim(5, 100)
        ax.set_title('Plasmon number vs. L, $\lambda = $' + lam)

        pdf.savefig(fig)
            
    pdf.close()

def ccplot():

    def plotcij():

        amax = 1.0

        #fig, ax = plt.subplots()
        figcij.subplots_adjust(right=0.87)
        cs = axcij[row][col].imshow(cij, cmap='bwr', vmin = -amax, vmax= amax)
        axcij[row][col].invert_xaxis()
        cax = figcij.add_axes([0.9, 0.1, 0.05, 0.8])
        figcij.colorbar(cs, cax=cax)

    def plotamn():

        amax = 1.0
        amn = phi.dot(cij).dot(phi.transpose())
        figamn.subplots_adjust(right=0.87)
        #fig, ax = plt.subplots()
        cs = axamn[row][col].imshow(amn, cmap='bwr', vmin = -amax, vmax= amax)
        axamn[row][col].invert_xaxis()
        cax = figamn.add_axes([0.9, 0.1, 0.05, 0.8])
        figamn.colorbar(cs, cax=cax)
        

    def finishcij():
        
        #plt.gca().invert_yaxis()
        
        figcij.suptitle(r'$ \langle {{{}}} | c^{{\dagger}}_i c_j | {{{}}} \rangle $ at $\lambda = {{{}}}$'.format( cur[0], cur[1], lam), size=50)
        pdfcij.savefig(figcij)
        print('closecij')
        #pdfcij.close()

    def finishamn():

        #plt.gca().invert_yaxis()
        figamn.suptitle(r'$ \langle {{{}}} | a^{{\dagger}}_m a_n | {{{}}} \rangle $ at $\lambda = {{{}}}$'.format( cur[0], cur[1], lam), size=50)
        pdfamn.savefig(figamn)
        print('closeamn')
        #pdfamn.close()


    def gen_phi():

        return np.array([[ np.sqrt(2/(L + 1)) * np.sin(n * np.pi * x/ (L + 1)) for x in range(1, L + 1)] for n in range(1, L + 1)])


    lam = sys.argv[1]


    prefix = syspre() + '/Desktop/Code/iTensor/collect/plasmon/high' + lam + '/cc/'
    dirs = os.listdir( prefix)
    dirs = sorted(dirs, key= lambda x: [ int(x.split('_')[-2]), int(x.split('_')[-1]), int(x.split('_')[1])])

    lengths = defaultdict(int)
    for dir in dirs:

        pair = ( int(dir.split('_')[-2]), int(dir.split('_')[-1]) )
        lengths[pair] += 1

    print(dirs)
    plotloc = syspre() + '/Desktop/Code/iTensor/plots/cc/'

    # init
    cur = ex = (int(dirs[0].split('_')[-2]), int(dirs[0].split('_')[-1]))
    cnt = 0
    dim = int(np.ceil(np.sqrt( lengths[ex])))

    figcij, axcij = plt.subplots(dim, dim, figsize =(dim *5, dim * 5))
    figamn, axamn = plt.subplots(dim, dim, figsize =(dim *5, dim * 5))

    pdfcij = matplotlib.backends.backend_pdf.PdfPages( plotloc + 'out_cij_lam{}.pdf'.format(lam))
    pdfamn = matplotlib.backends.backend_pdf.PdfPages( plotloc + 'out_amn_lam{}.pdf'.format(lam))

    for file in dirs:

        ex = (int(file.split('_')[-2]), int(file.split('_')[-1]))
        L = int(file.split('_')[1])

        print(ex, L)
        if ex != cur:

            finishcij()
            finishamn()
            #pdfcij = matplotlib.backends.backend_pdf.PdfPages( plotloc + 'out_cij{}_{}_{}.pdf'.format(lam, ex[0], ex[1]))
            #pdfamn = matplotlib.backends.backend_pdf.PdfPages( plotloc + 'out_amn{}_{}_{}.pdf'.format(lam, ex[0], ex[1]))
            cur = ex
            cnt = 0
            dim = int(np.ceil(np.sqrt( lengths[ex])))
            figcij, axcij = plt.subplots(dim, dim, figsize =(dim *5, dim * 5))
            figamn, axamn = plt.subplots(dim, dim, figsize =(dim *5, dim * 5))

        row = cnt // dim
        col = cnt % dim

        phi = gen_phi()
        cij = np.loadtxt( prefix + file)
        plotcij()
        plotamn()

        cnt += 1

    pdfcij.close()
    pdfamn.close()



def energytrendplot():

    def line(x, a, b):
        return a * x + b

    pdf = matplotlib.backends.backend_pdf.PdfPages( syspre() + '/Desktop/Code/iTensor/plots/energytrendplot.pdf' )
    plasmondir = syspre () + '/Desktop/Code/iTensor/collect/plasmon/'

    figgs, axgs = plt.subplots( figsize=(7, 8))
    fig1, ax1 = plt.subplots( figsize=(7, 8))
    figpl, axpl = plt.subplots( figsize =(7,8))
    axes = [axgs, ax1, axpl]
    key = ['GS', '1st ex', 'Plasmon']
    styles = ['dashdot', 'solid', 'dashed']
    gsfits=[]
    exfits = []
    lams = []

    for i, dir in enumerate(sorted( glob.glob( plasmondir  + 'high*'), key= lambda x: float(x[ len(plasmondir) + 4:]))):

        lam = dir[ len(plasmondir) + 4:]
        header = np.loadtxt(  dir + '/header')
        raw = np.loadtxt( dir + '/energy')

        refheader = np.log10(header)

        axgs.plot( refheader, - np.log10( -raw[:, 0]), label='$\lambda = $' + str(lam), linestyle=styles[i//10])
        ax1.plot( refheader, - np.log10( -raw[:, 1]), label='$\lambda = $'+ str(lam), linestyle=styles[i//10])
        axpl.plot( refheader, np.log10(raw[:, 1] - raw[:, 0]), label='$\lambda = $'+ str(lam), linestyle=styles[i//10])

        gsfit, cov = curve_fit(line, refheader,  np.log10( -raw[:, 0]))
        exfit, cov = curve_fit(line, refheader,  np.log10( -raw[:, 1]))

        gsfits += [gsfit[0]]
        exfits += [exfit[0]]
        lams += [ float(lam)]


    for i, ax in enumerate(axes):
        #ax.legend()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.4))

        if i < 2:
            ax.set_ylabel( r'$ - log_{{10}}(|E_{{{}}}|)$'.format(key[i]))

        else:
            ax.set_ylabel(r'$  log_{{10}}(|E_{{{}}}|)$'.format(key[i]))
        
        ax.set_xlabel(r'$log_{{10}} L$')
        ax.set_title(key[i] + ' Energy vs. L')
        #ax.set_xscale('log')
        #ax.set_xlim(5, 100)

    figfit, axfit = plt.subplots()

    axfit.plot(lams, gsfits, label='GS fit')
    axfit.plot(lams, exfits, label='1st ex fit')
    axfit.set_xlabel('$\lambda$')
    axfit.set_ylabel('Power law exponent')
    axfit.set_title('Power law exponent $(E \sim L^a)$ vs. $\lambda$')
    axfit.legend()

    figgs.tight_layout()
    fig1.tight_layout()
    figpl.tight_layout()
    figfit.tight_layout()

    pdf.savefig(figgs) 
    pdf.savefig(fig1)
    pdf.savefig(figpl)
    pdf.savefig(figfit)
    pdf.close()

def tcdplot():

    workdir = os.getcwd() + '/collect/plasmon/tcd/'
    plotdir = os.getcwd() + '/plots/'
    dirs = [0.5, 1.0]
    L = 98
    dis = 0
    upper = 21

    gap = 0.15
    for d in dirs:
        fig, ax = plt.subplots( figsize = (100, 15))
        for left in range(1, upper - 1):
            for right in range(left + 1, upper):

                tcd = np.loadtxt( workdir + '{}/TCD_{}_ex_{}_{}'.format(d, L, left, right))

                ref = np.arange(  1, L + 1) + (L + 5) * (left - 1)
                ax.scatter(ref, tcd + gap * (right - 2), label='GS to {}ex'.format(right - 1), s=1)
                ax.plot(ref, tcd + gap * (right - 2))

            ax.annotate( 'Transition to excited state {}'.format(left ) , ( -L , gap * ( left - 1)), size=20)
            ax.annotate( 'Transition from excited state {}'.format(left - 1), ( (L + 5) * (left - 1) + 2 , gap * ( left - 2)), size=20)
            
        ax.set_xlim( - L, L * (upper - 1))
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.4))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_title( 'TCD: GS to {} lowest Ex, L = {}, disorder = {}, $\lambda = ${}'.format(upper - 1, L, dis, d))
        fig.tight_layout()
        fig.savefig(plotdir + 'TCDGS{}_{}.pdf'.format(L, d))
        plt.show()


def gpiplot():

    def cal_GPI(tcd, ees):

        gpi = tcd.dot(ees).dot(tcd)

        #print(gpi)
        return gpi

    def ee(L, Lx, disx, disy, lam):
        
        arr = np.zeros((L, L))

        int_ee = 2 * lam
        zeta = 0.5
        z = 1.0
        ex = 0.2
        dis = [disx, disy]

        def distance():
            return np.sqrt( ( - dis[0][i] + xd + dis[0][j]) ** 2 + ( - dis[1][i] + yd + dis[1][j] ) ** 2)

        for i in range(L):
            for j in range(L):
                #ex = 0
                #zeta = 0
                xi = i % Lx
                yi = i // Lx 

                xj = j % Lx
                yj = j // Lx

                xd = xj - xi
                yd = yj - yi

                
                if abs(xd) + abs(yd) == 1:
                    factor = 1 - ex
                
                else:
                    factor = 1

                arr[i, j] =  z * int_ee * factor / ( distance() + zeta)

        return arr
        

    workdir = os.getcwd() + '/collect/plasmon/tcd/'
    plotdir = os.getcwd() + '/plots/'
    dirs = [1.0, 0.5]
    L = 98
    dis = 0
    upper = 21

    disx = np.zeros(L)
    disy = np.zeros(L)
    
    gpis = np.zeros((upper - 2, upper - 2))
    for d in dirs:
        ees = ee(L, L, disx, disy, d)
        fig, ax = plt.subplots( figsize = (15, 15))
        for left in range(1, upper - 1):
            for right in range(left + 1, upper):

                tcd = np.loadtxt( workdir + '{}/TCD_{}_ex_{}_{}'.format(d, L, left, right))
                gpi = cal_GPI(tcd, ees)
                
                gpis[right - 2][left - 1] = gpi

        print(gpis)
        cs = ax.imshow(gpis, origin='lower', extent=[0, upper - 2, 1, upper - 1])
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title( 'GPI: GS to {} lowest Ex, L = {}, disorder = {}, $\lambda = ${}'.format(upper - 1, L, dis, d))
        ax.set_xlabel('From state (0 = GS)')
        ax.set_ylabel('To state')
        fig.colorbar(cs, ax=ax)
        fig.tight_layout()
        fig.savefig(plotdir + 'GPIGS{}_{}.pdf'.format(L, d))
        plt.show()

        
if __name__ == '__main__':
    
    lengthcompplot()
    #weightplot()
    #longplot()
    #plot2d()
    #plasmonplot()
    #highplasmonplot()
    #ccplot()
    #energytrendplot()
    #extrendplot()
    #tcdplot()
    #gpiplot()