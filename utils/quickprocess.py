import glob
import numpy as np


vals = glob.glob('log_new_krylov/energy*')


for val in vals:
    
    cur = val[ len('log_new_krylov/energy'):]

    dim, up, dn, t = cur.split('_')
    print(dim, up, dn, t)

    temp = np.loadtxt(val)
    np.savetxt( 'log_new_krylov/energy{}x{}_{}_{}_{}'.format(dim, dim, t, up, dn), temp)