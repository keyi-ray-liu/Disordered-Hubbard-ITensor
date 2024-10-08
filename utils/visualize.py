import matplotlib.pyplot as plt
import networkx as nx
import sys
from matplotlib import rc
import numpy as np
import glob
import os
from collections import defaultdict

def format_pair(vals:list[str]):

    pairs = []
    pairs.append(float(vals[0]))

    for val in vals[1:]:
        val = val.split(',')[0]
        op, num = val.split('(')
        num = int(num)

        pairs.append([op, num])
    return pairs

def work(system ="DPT", dir=os.getcwd(), args=dict(), plot=True):

        
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 


    #G = nx.Graph()

    #dir = sys.argv[1]

    if not args:

        if system == "DPT":
            files = glob.glob( dir + '/HS*')
            LR = np.loadtxt( dir + '/LR')
            col = 3

            if 'includeUFalse' in dir:

                mid = LR.shape[0] // 2
                LR = np.concatenate((LR[:mid] , np.zeros(4) , LR[mid:]))

        elif system == "QE" :
            files = glob.glob( dir + '/H*')
            col = 3

        elif system == "SD" :
            files = glob.glob( dir + '/H*')
            col = 2


        out = 'plots/' 

        fig, axes = plt.subplots(len(files), col, figsize=(10 * col, 5 * len(files)))

        if len(files) == 1:
            axes = [axes]

    else:

        files = args["files"]
        fig = args["fig"]
        axes = args["axes"]

    if not plot:
        minval = 10000


    for i, file in enumerate(files):

        
        with open(file, 'r') as f:
            raws = f.readlines()

    
        rawpairs = [ raw.strip().split(' ') for raw in raws[1:-1]]

        if system == "DPT":

            total = LR.shape[0] + 2

            nonU = np.zeros((total, total))
            U = np.zeros((total, total))

            L = []
            R = []
            lower = []
            upper = [total, 0.0]
            spatial = []

        elif system == "QE":

            potential = []
            N = np.loadtxt(dir + '/occ').shape[-1]
            ee = np.zeros((N, N))
            qe = np.zeros((N, N))

        elif system == "SD":

            potential = []
            try:
                N = np.loadtxt( dir + '/occ').shape[-1]

            except OSError as e:
                N = np.loadtxt( dir + '/occup').shape[-1]

            ee = np.zeros((N, N))
            qe = np.zeros((N, N))

        print(rawpairs)
        for raw in rawpairs:
            
            pair = format_pair(raw)
            # these are the onsite terms
            if len(pair) == 2 and pair[1][0] != "I":
                coupling, ops= pair
                _, s = ops

                if system == "DPT":
                    if s == total - 1:
                        lower = [s, coupling]

                    elif s == total:
                        upper = [s, coupling]

                    elif LR[s - 1] > 0:
                        L.append([s, coupling])

                    elif LR[s - 1] < 0:
                        R.append([s, coupling])

                    else:
                        spatial.append([s, coupling])

                elif system == "QE":

                    potential.append([s, coupling])

                elif system == "SD":

                    potential.append([s, coupling])


            elif len(pair) == 3:
                coupling, op1, op2= pair

                _, s1 = op1
                _, s2 = op2

                if system == "DPT":
                    nonU[ s1 -1, s2- 1] += coupling

                elif system == "QE":
                    ee[ s1 - 1, s2 - 1] += coupling
                    
                    if not plot and  s1 < 70 and s2 > 90:
                        minval = min(minval, s1, s2)


                elif system == "SD":
                    ee[ s1 - 1, s2 - 1] += coupling
                    


            elif len(pair) == 4:
                coupling, op1, op2, op3= pair

                _, s1 = op1
                _, s2 = op2
                _, s3 = op3

                #print(coupling)

                if system == "DPT":
                    U[ s2 -1, s3- 1] += coupling

                elif system == "QE":
                    qe[s1 -1, s3 -1] += coupling

        if not plot:
            print(minval)

                
        if plot:
            s = 1
            ax : plt.Axes = axes[i][0]

            if system == "DPT":
                ax.scatter( np.array(L)[:, 0], np.array(L)[:, 1], c='red', label='L', s=s)
                ax.scatter( np.array(R)[:, 0], np.array(R)[:, 1], c='blue', label='R', s=s)

                if spatial:
                    ax.scatter( np.array(spatial)[:, 0], np.array(spatial)[:, 1], c ='black', label='spatial', s=s)
                ax.scatter( lower[0], lower[1], c='orange', label='lower', s=s)
                ax.scatter( upper[0], upper[1], c='magenta', label='upper', s=s)
                ax.set_ylim(-5, 5)

            elif system == "QE":
                ax.scatter( np.array(potential)[:, 0], np.array(potential)[:, 1], c='red', label='L', s=s)

            elif system == "SD":
                ax.scatter( np.array(potential)[:, 0], np.array(potential)[:, 1], c='red', label='L', s=s)


                #     G.add_node( s, onsite=coupling)

                # elif len(pair) == 3:

                #     coupling, ops1, ops2, = pair
                #     op1, s1 = ops1
                #     op2, s2 = ops2

                #     G.add_edge( s1, s2, coupling=coupling)

            ax.legend()
            ax.set_xlabel('Site numer')
            ax.set_ylabel('Onsite potential')
            ax.set_title('Onsite interactions vs. Site, ' + file)
            
            #nx.draw_circular(G, with_labels=True, font_weight='bold')


            ax : plt.Axes = axes[i][1]

            if system == "DPT":
                im = ax.imshow( nonU, interpolation='none', rasterized=True, vmax=2.0, vmin=-2.0)
                ax.set_title(r'$c_k c_l $ two-body interaction')

            elif system == "QE" :
                im = ax.imshow( ee, interpolation='none', rasterized=True,)
                ax.set_title(r'$e-e$ and hopping interaction for site pair $i, j$')

            elif system == "SD" :
                im = ax.imshow( ee, interpolation='none', rasterized=True,
                               #norm='symlog', 
                               extent=[1, ee.shape[0], 1, ee.shape[0] ])
                ax.set_title(r'$e-e$ and hopping interaction for site pair $i, j$')

            ax.invert_yaxis()
            fig.colorbar(im)

            if system != 'SD':
                ax : plt.Axes = axes[i][2]

                if system == "DPT":
                    im2 = ax.imshow( U, interpolation='none', rasterized=True)
                    ax.set_title(r'$ N_{lower} c_k c_l $ three-body interaction')

                elif system == "QE":
                    im2 = ax.imshow( qe, interpolation='none', rasterized=True,)
                    ax.set_title('QE interaction')

                ax.invert_yaxis()
                fig.colorbar(im2)


    if not args and plot:
        fig.tight_layout() 
        fig.savefig(out + 'visualize{}.pdf'.format(dir), dpi=1000)



def graph():

    def get_pos():

        row = col = 3
        mag = 2
        pos = {1: (mag, 0)}
        for i in range(row):
            for j in range(col):

                ii = i * mag
                jj = j * mag
                pos[ i * row + j + 2] = (ii, jj + mag)
        pos[11] = (mag, ( row + 1) * mag )
        return pos

    file = sys.argv[1]

    with open(file, 'r') as f:
        raws = f.readlines()
    rawpairs = [ raw.strip().split(' ') for raw in raws[1:-1]]

    G = nx.Graph()
    hop = defaultdict(float)
    den = defaultdict(float)
    onsite = defaultdict(float)

    for raw in rawpairs:
        
        pair = format_pair(raw)
        # these are the onsite terms
        if len(pair) == 2 and pair[1][0] != "I":
            coupling, ops= pair
        
            _, s = ops

            if coupling != 0:
                onsite[ s] += coupling

        elif len(pair) == 3:
            coupling, op1, op2= pair
            id1, s1 = op1
            id2, s2 = op2
            s1, s2 = sorted((s1, s2))
            
            if 'C' in id1 or 'C' in id2:
                if coupling != 0:
                    hop[ (s1 , s2  )] += coupling

            else:

                if coupling != 0:
                    den[ (s1, s2) ] += coupling

        else:
            raise ValueError("We do not support more than 2 op at the moment")
        
    hop = { val : hop[val]/4 for val in hop}
    
    for val in onsite:
        G.add_node( val)

    for val in hop:
        G.add_edge( *val)
    
    for val in den:
        G.add_edge( *val)

    pos = get_pos()

    plt.figure()
    nx.draw(
    G, pos, edge_color='black', width=1, linewidths=1,
    node_size=500, node_color='pink', alpha=0.9,
    labels={ val:  '{:.3g}'.format(onsite[val])  for val in onsite},
    connectionstyle='arc3,rad=0.2', arrows=True
    )

    nx.draw_networkx_edge_labels(
    G, pos,
    edge_labels={val :  '{:.3g}'.format(hop[val])  for val in hop},
    font_color='red',
    connectionstyle='arc3,rad=0.2'
    )

    nx.draw_networkx_edge_labels(
    G, pos,
    edge_labels={val : '{:.3g}'.format(den[val]) for val in den},
    font_color='orange',
    connectionstyle='arc3,rad=0.2'
    )
    plt.show()

if __name__ == '__main__':

    work( system = sys.argv[2], dir=sys.argv[1])

    #graph()

