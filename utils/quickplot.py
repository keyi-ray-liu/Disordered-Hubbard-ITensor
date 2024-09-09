import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
from matplotlib import rc
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.append(os.path.abspath("/Users/knl20/Desktop/Code/TN/LT"))
from searchkey import *
from matplotlib.gridspec import GridSpec
from numpy.fft import rfft
from scipy.signal import find_peaks
from matplotlib.backends.backend_pdf import PdfPages


def load_data(raw, tswitch):
    raw_data = np.loadtxt(raw)

    time = raw_data[:, 0]
    occ = raw_data[:, 5]
    cur = raw_data[:, 3]

    idstop = np.argwhere(time <= tswitch * 2).flatten()[-1] + 1
    time = time[:idstop]
    occ = occ[:idstop]
    cur = cur[:idstop]

    return time, occ, cur


def sign_change(arr, avg_step=200):

    signs = np.sign(arr)
    signs = np.roll(signs, 1) - signs
    
    idx = np.argwhere(signs != 0).flatten()[-2:]

    if len(idx) < 2:
        idx = [-avg_step, -1]
    return idx


def MF1(file):

    raw = np.loadtxt(file)

    time = np.concatenate(([0.0], raw[:,0]))

    #cur1 = np.concatenate(([0.0] , raw[:,2]))
    cur2 = np.concatenate(([0.0] , raw[:,3]))
    #cur3 = np.concatenate(([0.0] , raw[:,4]))

    occ = np.concatenate(([1.0] ,raw[:,5]))

    #avg = np.mean( raw[:, 1:4], axis=1)

    fig ,ax = plt.subplots(figsize=(15, 5))

    #L = searchkey('N', file)
    #U = searchkey('U', file)


    ax : plt.Axes 
    #ax.plot(time, cur1, label='cur1')
    #ax.plot( time, (cur3 - cur1)/ ( 2 *np.pi), label='cur3 - cur1 ')
    #ax.plot( time, (cur3 + cur1 + cur2)/3, label='sum ')
    ax.plot(time, cur2, label='cur2')
    #ax.plot(time, cur3, label='cur3')
    #ax.plot(time, avg, label='avg')

    ax2 : plt.Axes = ax.twinx()
    ax2.scatter(time, occ, s=1, label='dd', c='orange')

    #color = 'red'
    #ref = np.loadtxt('/Users/knl20/Desktop/Code/TN/non-interacting/{}currentCC{}'.format(L, U))

    #reftime = np.arange(0, 32, 0.25)

    # ax.plot( reftime, ref[:reftime.shape[0]] , label=r'$I_{L\rightarrow R}$, Stage 2, non-interacting exact', color=color,
    #         linestyle='dotted'
    #         )
    ax.legend()
    ax2.legend()
    ax2.set_ylim(0, 1.1)

    ax.set_ylabel('Current')
    ax2.set_ylabel('DD ')

    fig.savefig( '../plots/MFtest{}.pdf'.format(file[:-len('.txt')]))
    return fig


def test():

    def format_axes(axs):
        for k in axs:
            axs[k].text(0.5, 0.5, f"ax: {k}", va="center", ha="center")
            axs[k].tick_params(labelbottom=False, labelleft=False)


    fig, axs = plt.subplot_mosaic([["phocc", "occ1", "cur1", "phcur"], ["phocc", "occ2", "cur2", "phcur"]], constrained_layout=True,
                                gridspec_kw={'width_ratios':[1, 1, 1, 1],
                                'height_ratios':[1, 1]})
    
    print(axs)
    fig.suptitle("subplot_mosaic")
    format_axes(axs)
    plt.show()



def diff(occ):

    #ori = np.linspace(0,5, 100)
    #gap = 5/100

    f = np.gradient(occ)
    plt.plot(f)

    plt.show()


def test_fft():

    t = 0.3
    x = np.arange(500) * t
    y = np.sin( 0.1 * 2 * np.pi * x) + 0.5 * np.sin( 0.2 * 2 * np.pi * x) 

    xf = np.abs(rfft(y))

    print(np.argsort(xf)[::-1][:10])
    plt.plot(np.arange(1, xf.shape[0]), xf[1:])
    plt.show()

def get_fft(occ, sample_time):

    tstep = sample_time[1] - sample_time[0]

    print("tstep: {}".format(tstep))
    y = occ
    yf = rfft(y)

    yf  = np.abs(yf)

    peaks, _ = find_peaks(yf[1:], prominence=0.14 * (np.amax(yf[1:]) - np.amin(yf[1:]) ))
    print(peaks)

    #grad = np.gradient(yf[1:])
    #diff_grad = grad[1:] - grad[:-1]
    #idx_grad = np.argsort(diff_grad)[:2]
    
    xf = np.arange(yf.shape[0]) / occ.shape[0] / tstep

    # all = occ.shape[0]
    # occ = occ[all//2:]
    # N = occ.shape[0]//2
    # amp = rfft(occ)

    # xf = np.linspace(0, N//2, N//2)
    # yf = 2.0/N * np.abs(amp[:N//2])



    return xf[1:], yf[1:], peaks

def phase_diagram(vs =0.125, fftseparate=True, fft=True):


    def plot_occ():
        
        time, occ, _ = load_data(raw, tswitch)

        # [left, bottom, width, height
        axocc :plt.Axes = axes[ 'occ' + str(occ_cnt)]
        con = ConnectionPatch(xyA=(0.05, 0.5), xyB=(float(mu), float(U)), coordsA="axes fraction", coordsB="data", axesA=axocc, axesB=axphocc, color="black")
        con.set_in_layout(False)
        axocc.add_artist(con)
        
        #ref = np.arange(occ.shape[0])/step
        #axocc.tick_params(left=False, right=True, labelleft=False, labelright=True)

        #vy = occ
        #vx = ref

        axocc.plot(time, occ, color=plotted_occ[key])
        #axocc.plot( time, grad, color=plotted_occ[key])

        axocc.set_ylim(0, 1)
        lo, hi = axocc.get_ylim()
        l, r = axocc.get_xlim()
        
        axocc.vlines( [ t], [lo], [hi], linestyles='dotted')
        axocc.hlines( [ vs], [l], [r], linestyles='dashed')
    
        axocc.set_title(r'$\langle n_1 \rangle: \mu$' + ' = {}, U = {}'.format(mu, U), fontsize=fontsize)
        #axocc.legend()

        axocc.set_xlabel('t $(1/\omega_0)$', fontsize=fontsize)
        axocc.set_ylabel(r'$\langle n_1\rangle$', fontsize=fontsize)
        axphocc.scatter([float(mu)], [float(U)], c='black', s=s)

        axocc.tick_params(axis='both', labelsize=tickfontsize)

    def fft_occ():

        time, occ, _ = load_data(raw, tswitch)

        idswitch = np.argwhere( time <= tswitch ).flatten()[-1] + 1

        # [left, bottom, width, height

        if fftseparate:
            axfft :plt.Axes = axes[ 'occ' + str(occ_cnt)]
            con = ConnectionPatch(xyA=(0.05, 0.5), xyB=(float(mu), float(U)), coordsA="axes fraction", coordsB="data", axesA=axfft, axesB=axphocc, color="black")
            con.set_in_layout(False)
            axfft.add_artist(con)

        else:
            axfft : plt.Axes = axes[ 'offt' + str(occ_cnt)]



        vx, vy, idxpeak = get_fft(occ[idswitch:], time[idswitch:])
        axfft.plot(vx, vy, color=plotted_occ[key], label='S3')
        
        for peak in idxpeak:
            axfft.scatter( vx[peak], vy[peak], label='{:.3g}'.format(vx[peak]))
    
        axfft.set_title(r'$f_{{FFT}}(\langle n_1 \rangle) : \mu$' + ' = {}, U = {}'.format(mu, U), fontsize=fontsize)
        #axocc.legend()

        axfft.legend(fontsize=fontsize)
        axfft.set_xlim(0, axfft.get_xlim()[1]/4)
        #axfft.set_xlabel(r'$ \omega = \frac{{k}}{{N\Delta t}} $', fontsize=fontsize)
        axfft.set_xlabel('$ \omega / (N\Delta t) $', fontsize=fontsize)
        axfft.set_ylabel(r'$|\sum_{m=0}^{N-1} \langle n_1\rangle e^{-2\pi i\omega m/N } |$', fontsize=fontsize)
        axphocc.scatter([float(mu)], [float(U)], c='black', s=s)


        axfft.tick_params(axis='both', labelsize=tickfontsize)

    def plot_cur():

        time, _, cur = load_data(raw, tswitch)

        axcur :plt.Axes = axes['cur' + str(cur_cnt)]

        con = ConnectionPatch(xyA=(0.95, 0.5), xyB=(float(mu), float(U)), coordsA="axes fraction", coordsB="data", axesA=axcur, axesB=axphcur, color="black")
        con.set_in_layout(False)
        axcur.add_artist(con)

        #axcur.tick_params(left=True, right=False, labelleft=True, labelright=False)
        axcur.plot(time, cur, color=plotted_current[key])
        axcur.vlines( [ t], [axcur.get_ylim()[0]], [axcur.get_ylim()[1]], linestyles='dotted')



        axcur.set_title(r'$I/\mu$: $\mu$' + '= {}, U = {}'.format(mu, U), fontsize=fontsize)
        #ax.legend()
        axphcur.scatter([float(mu)], [float(U)], c='black', s=s)

        axcur.set_xlabel('t $(1/\omega_0)$', fontsize=fontsize)
        axcur.set_ylabel('$I/\mu$', fontsize=fontsize)
        
        axcur.tick_params(axis='both', labelsize=tickfontsize)
    

    def fft_cur():

        time, _, cur = load_data(raw, tswitch)

        idswitch = np.argwhere( time <= tswitch ).flatten()[-1] + 1

        # [left, bottom, width, height

        if fftseparate:
            axfft :plt.Axes = axes[ 'cur' + str(cur_cnt)]
            con = ConnectionPatch(xyA=(0.95, 0.5), xyB=(float(mu), float(U)), coordsA="axes fraction", coordsB="data", axesA=axfft, axesB=axphcur, color="black")
            con.set_in_layout(False)
            axfft.add_artist(con)

        else:
            axfft :plt.Axes = axes[ 'cfft' + str(cur_cnt)]

        vx, vy, idxpeak = get_fft(cur[idswitch//2:idswitch], time[idswitch//2:idswitch])
        axfft.plot(vx, vy, color=plotted_current[key], linestyle='dotted', label='Stage 2 ($t>t_s/2$)')

        for peak in idxpeak:
            axfft.scatter( vx[peak], vy[peak], label='S2: {:.3g}'.format(vx[peak]) )
        
        vx, vy, idxpeak= get_fft(cur[idswitch:], time[idswitch:])
        axfft.plot(vx, vy, color=plotted_current[key], label='Stage 3')

        for peak in idxpeak:
            axfft.scatter( vx[peak], vy[peak], label='S3: {:.3g}'.format(vx[peak]), marker='x')

        axfft.legend(fontsize=fontsize)
    
        axfft.set_title(r'$I/\mu : \mu$' + ' = {}, U = {}'.format(mu, U), fontsize=fontsize)
        #axocc.legend()

        axfft.set_xlim(0, axfft.get_xlim()[1]/4)
        #axfft.set_xlabel(r'$ \omega = \frac{{k}}{{N\Delta t}} $', fontsize=fontsize)
        axfft.set_xlabel('$ \omega/ (N\Delta t) $', fontsize=fontsize)
        axfft.set_ylabel(r'$|\sum_{m=0}^{N-1} I e^{-2\pi i\omega m/N } |$', fontsize=fontsize)
        axphocc.scatter([float(mu)], [float(U)], c='black', s=s)


        axfft.tick_params(axis='both', labelsize=tickfontsize)

    modes =[ "discrete", "equal"]



    for mode in modes:


        if fftseparate:
            file = '../plots/MFphase{}mode{}fft{}.pdf'.format(vs, mode, fft)
            mosaic = [
                        ["phocc", "occ1", "cur1", "phcur"], 
                        ["phocc", "occ2", "cur2", "phcur"], 
                        ["phocc", "occ3", "cur3", "phcur"]
                                                ]
        else:
            mosaic = [
                        ["phocc", "occ1", "offt1", "cfft1", "cur1", "phcur"], 
                        ["phocc", "occ2", "offt2", "cfft2", "cur2", "phcur"], 
                        ["phocc", "occ3", "offt3", "cfft3", "cur3", "phcur"]
                                                ]
            file = '../plots/MFphase{}mode{}.pdf'.format(vs, mode)
    
        with PdfPages(file)  as pdf:
            factor_str = {1 : '$U, U, U, U$',
                        2 : r'$\frac{U}{2}  , \ U,  U,   \ \frac{U}{2}$'}
                        

            plotted_occ = {
                            ('5.0', '1.0') if vs == 1.0 else ('2.0', '1.0'): 'red',
                            ('2.0', '1.0') if vs == 1.0 else ('1.0', '1.0'): 'blue',
                            ('1.0', '1.0') if vs == 1.0 else ('0.5', '1.0'): 'black',
                            }
            
            plotted_current = {
                            ('5.0', '1.0') if vs == 1.0 else ('2.0', '1.0'): 'red',
                            ('2.0', '1.0') if vs == 1.0 else ('1.0', '1.0'): 'blue',
                            ('1.0', '1.0') if vs == 1.0 else ('0.5', '1.0'): 'black',
                            }
            
            #eachrow = max( len(plotted_current.keys()), len(plotted_occ.keys()))

            
            for factor in [1, 2]:
                
                tickfontsize = 20
                fontsize=20
                s = 5
                avg_step = 200
                Ns = [ 258]

                if vs == 1.0:
                    ts = np.array([1/4])

                elif vs == 0.125:
                    ts = np.array([1/2])

                # fig = plt.figure(
                #     figsize = (4 * 4, 6 * eachrow),
                #     layout="constrained"
                #     )
                # gs = GridSpec(eachrow, 4, figure=fig) 

                fig, axes = plt.subplot_mosaic(mosaic, constrained_layout=True,
                                    gridspec_kw={'width_ratios':[1 for _ in range(len(mosaic[0]))],
                                    'height_ratios':[0.7 for _ in range(len(mosaic))]}, figsize=(6 * len(mosaic[0]), 16))


                #dirs = "examples_v{}/".format(vs)
                dirs = "new/"

                for i, n in enumerate(Ns):
                    
                    tswitch = ts * n - 1
                    for j, t in enumerate(tswitch):
                        
                        #axes[i][j * metric + 1].set_axis_off()

                        raws = glob.glob( dirs + '*N{}*tswitch{}*factor{}*{}*'.format( n, float(t), factor, mode))


                        raws = sorted(raws, reverse=True)

                        mus = sorted(list(set([searchkey('mu', f) for f in raws])))
                        mus = { mu : k for k, mu in enumerate(mus)}

                        Us = sorted(list(set([searchkey('U', f) for f in raws])))
                        Us = { U : k for k, U in enumerate(Us)}

                        occ_data = np.zeros((len(mus), len(Us)))
                        cur_data = np.zeros((len(mus), len(Us)))
                        #diffs = np.zeros((len(mus), len(Us)))

                        extent = [ float(min(mus.keys())), float(max(mus.keys())), float(max(Us.keys())), float(min(Us.keys()))]


                        occ_cnt = 1
                        cur_cnt = 1


                        for raw in raws:
                            
                            mu = searchkey( 'mu', raw)
                            U = searchkey('U', raw)
                        
                            _, occ, cur = load_data(raw, tswitch)
                            cur = cur / float(mu)

                            idx1, idx2 = sign_change(np.gradient(occ))

                            #print(idx1, idx2)
                            occ_to_avg = occ[idx1:idx2]
                            cur_to_avg = cur[-avg_step:]

                            #diff = np.amax(occ_to_avg) - np.amin(occ_to_avg)
                            occ_avg = np.average(occ_to_avg)
                            cur_avg = np.average(cur_to_avg)

                            occ_data[ mus[mu], Us[U]] = occ_avg
                            cur_data[ mus[mu], Us[U]] = cur_avg



                        axphocc : plt.Axes = axes['phocc']
                        im = axphocc.imshow(occ_data.transpose(), extent=extent, cmap='bwr', vmin=0, vmax=1)
                        axphocc.invert_yaxis()
                        axphocc.set_xlabel('$\mu$', fontsize=fontsize)
                        axphocc.set_ylabel('U', fontsize=fontsize)
                        axphocc.set_title(r'$\overline{\langle n_1\rangle}$: ' + ' N = {}, tswitch = {}'.format(n, t), fontsize=fontsize)

                        divider = make_axes_locatable(axphocc)
                        cax = divider.append_axes("bottom", size="3%", pad=0.65)

                        cb = fig.colorbar(im, cax=cax, orientation="horizontal"
                                    )
                        cb.ax.tick_params(labelsize=tickfontsize)
                        cb.set_label(label=r'$\overline{\langle n_1\rangle}$', fontsize=fontsize)
    
                        
                        axphocc.tick_params(axis='both', labelsize=tickfontsize)

                        axphcur : plt.Axes = axes['phcur']
                        im = axphcur.imshow(cur_data.transpose(), extent=extent, cmap='autumn')
                        axphcur.invert_yaxis()
                        axphcur.set_xlabel('$\mu$', fontsize=fontsize)
                        axphcur.set_ylabel('U', fontsize=fontsize)
                        axphcur.set_title(r'$\overline{I/\mu}$: ' + ' N = {}, tswitch = {}'.format(n, t), fontsize=fontsize)

                        divider = make_axes_locatable(axphcur)
                        cax = divider.append_axes("bottom", size="3%", pad=0.65)

                        cb = fig.colorbar(im, cax=cax, orientation="horizontal"
                                    )
                        cb.ax.tick_params(labelsize=tickfontsize)
                        cb.set_label(label=r'$\overline{I/\mu}$', fontsize=fontsize)

                        axphcur.tick_params(axis='both', labelsize=tickfontsize#, left=False, right=True, labelleft=False, labelright=True
                                            )

                        
                        for raw in raws:
                            mu = searchkey( 'mu', raw)
                            U = searchkey('U', raw)
                            key = (U, mu)
                            
                            
                            if key in plotted_occ:
                                
                                if not fftseparate or not fft:
                                    plot_occ()
                                
                                if not fftseparate or fft:
                                    fft_occ()

                                occ_cnt += 1
                                

                            if key in plotted_current:
                                
                                if not fftseparate or not fft:
                                    plot_cur()

                                if not fftseparate or fft:
                                    fft_cur()

                                cur_cnt += 1
                            occ_data[ mus[mu], Us[U]]= occ_avg
                            cur_data[ mus[mu], Us[U]]= cur_avg 
                            #diffs[ mus[mu], Us[U]]= diff    





                fig.suptitle('$v_S ={}$, QPC: {} interaction'.format(vs, factor_str[factor]), fontsize=50)
                #fig.subplots_adjust(wspace=0, hspace=0)
                #fig.tight_layout()
                fig.savefig('factor{}.pdf'.format(factor))
                #plt.show()
                pdf.savefig(fig)
                    
                    


if __name__  == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 

    #MF1(glob.glob('*equal')[0])
    phase_diagram(fftseparate = False, fft=True)
    #test()

    # files = sys.argv[1:]

    # if not len(files):
    #     files = glob.glob('Current*')
    #     files = sorted( files, key=lambda x: [float(searchkey('U', x)), int(searchkey('N', x)), float(searchkey('tswitch', x))])

    # if len(files) > 1:

    #     with PdfPages('../plots/MFmulti.pdf') as pdf:
    #         for file in files:
    #             fig = MF1(file)
    #             fig.suptitle(file)

    #             pdf.savefig(fig)

    #test_fft()
    # else:

    #     fig = MF1(files)
    #     fig.suptitle(files)
    #     fig.savefig('../plots/MFtest'+ files)
        
