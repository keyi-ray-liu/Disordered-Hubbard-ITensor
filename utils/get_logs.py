import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import matplotlib.backends.backend_pdf
import os

def readlog(log_name, key, tag):

    logs = glob.glob('{}/res*'.format(log_name))
    for log in logs:

        print(log)
        with open(log, 'r') as f:
            raw = f.read()

        inds = [m.start() for m in re.finditer(key, raw)]
        #arr = [ float( raw[ind + len(key) : ind + len(key) + 12 ] ) for ind in inds]

        arr = [ float( raw[ind + len(key) : ind + raw[ind : ind + 30].find(' ') ] ) for ind in inds]

        arr = np.array(arr)
        newtag = log.replace('res', tag)
        np.savetxt( newtag, arr )





def prediction(log_name, key):

    vals = glob.glob('{}/{}*'.format(log_name, key))
    target_dir = os.getcwd() + '/predict/'
    for val in vals:
        
        trueval = val[ len(log_name) + len(key) + 1:]
        print(trueval)
        dim, t, up, dn =trueval.split('_')

        raw = np.loadtxt(val)

        gap = raw[-100:] - raw[-101:-1]
        np.savetxt( target_dir+ 'gap{}_{}_{}_{}'.format(dim, t, up, dn), gap)




def plot(log_name, key, exclude=[], begin=0):

    vals = glob.glob('{}/{}*'.format(log_name, key))
    cwd = os.getcwd() + '/'
    pdf = matplotlib.backends.backend_pdf.PdfPages( cwd + 'new_trend_{}_{}.pdf'.format(log_name, key) ) 
    uniq = np.infty
    

    for val in sorted(vals, key=lambda x: x.split('_')[1]):
        
        trueval = val[ len(log_name) + len(key) + 1:]

        print(trueval)
        dim, t, up, dn =trueval.split('_')
        cur = '_'.join([dim, up, dn])

        t = float(t)
        colors = {
            '4x4': 'blue', 
            '5x5':'red', 
            '6x6':'orange',
            '7x7': 'green',
            '3x3x3':'purple',

            'chain13': 'green',
            'chain31': 'red',
            'chain61': 'blue'
        }

        if cur not in exclude:

            if t != uniq:
                
                if uniq < np.infty:
                    ax.set_title( " |t/U| ={}".format( float(uniq) / 4))
                    ax.set_xlabel('Iteration number')
                    ax.set_ylabel(key)
                    ax.legend()
                    fig.tight_layout()
                    pdf.savefig(fig)
                    
                fig, ax = plt.subplots()
                uniq = t

            raw = np.loadtxt(val)
            raw = raw[begin:]
            ax.plot( np.arange(begin + 1, begin + len(raw) + 1), raw, label=cur, c=colors[dim], linewidth=0.02)

    ax.set_title( " |t/U| ={}".format( float(uniq) / 4))
    ax.set_xlabel('Iteration number')
    ax.set_ylabel('GS energy')
    ax.legend()
    fig.tight_layout()
    pdf.savefig(fig)
    pdf.close()

if __name__ == '__main__':

    exclude = []
    key = 'energy'
    log_name = 'log'

    readlog(log_name, 'maxlinkdim=', 'bond')
    readlog(log_name, 'energy=', 'energy')
    plot(log_name, key, exclude, 0)
    #prediction(log_name, key)


