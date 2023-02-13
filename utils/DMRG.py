import numpy as np
from numpy.linalg import svd
from scipy.sparse.linalg import eigsh
from copy import deepcopy


# defines the MPO operators on each site
# spin up and down

def MPOs():
    return [MPO(0), MPO(1), MPO(-1)]

def MPO(pos):

    t = -1.0
    
    if pos==1:
        # indices: phy L, phy R, b L, b R
        W = np.zeros((2, 2, 4, 4))

        # I
        for i in range(2):
            W[i, i, 0, 0] = 1
            W[i, i, -1, -1] = 1

        # c    <0|c|1> != 0
        W[0, 1, 1, 0] = 1
        # c+
        W[1, 0, 2, 0] = 1
        # tc
        W[0, 1, -1, 2] = t
        # tc+
        W[1, 0, -1, 1] = t

    else:
        W = np.zeros((2, 2, 4))
        # 0: start
        if pos == 0:

            for i in range(2):
                W[i, i, -1] = 1

            # tc
            W[0, 1, 2] = t
            # tc+
            W[1, 0, 1] = t

        # -1: end
        elif pos == -1:

            # I
            for i in range(2):
                W[i, i, 0] = 1

            # c+
            W[1, 0, 2] = 1
            # c
            W[0, 1, 1] = 1


    print(W)
    return W

# setting up the initial 'block'
def init_state(H):

    # random coeff
    psi = np.random.random((2, 2))
    U, D, V = svd(psi)


    Hleft, _, Hright = H
    # no transpose because U and V are technically 1D
    #left
    L = np.einsum('ij,ikl,km->jlm', U, Hleft, U)

    #right: indices are flipped for later consistency in SVD
    R = np.einsum('ji,ikl,mk->jlm', V, Hright, V)

    return L, R

# contract the current TN and solve the eigenproblem
def solve(L, R, H, phy_dim, max_bond):
    
    _, Hmid, _ = H

    #form the larger Hamiltonian
    # ijk: L: ( psi bond, H bond, psi bond)
    # abjm: Hmid dim: ( phy, phy, left H, right H)
    # cdmn: Hmid dim: ( phy, phy, left H, right H))
    # lnp: R: ( psi bond, H bond, psi bond)
    # Out: iaclkbdp: ( b, p, p, b (up) b, p, p, b (down))
    print('L dim: {}'.format(L.shape))
    print('R dim: {}'.format(R.shape))

    Hprime = np.einsum('ijk,abjm,cdmn,lnp->iaclkbdp',L,Hmid,Hmid,R)
    print('Hprime shape: {}'.format(Hprime.shape))

    a, b, c, d, _, _, _, _ = Hprime.shape

    if b != phy_dim or c != phy_dim:
        raise Exception('Error: physical dimension mismatch!')

    Hprime = Hprime.reshape(a * phy_dim * phy_dim * d, a* phy_dim * phy_dim * d, order='F')
    print('Hprime hermitian?', np.allclose(Hprime, Hprime.transpose()))

    # solve for the smallest eigenvalue
    w, v = eigsh(Hprime, k=1, which='SA')

    # now v is of dim b * b * d * d, we perform another SVD
    v = v.reshape(a * phy_dim, phy_dim * d, order='F')
    U, D, V = svd(v)

    print(D)

    print('U, V shape {}, {}'.format(U.shape, V.shape))

    # truncation
    if U.shape[1] > max_bond:
        V = V[:max_bond, :]
        U = U[:, :max_bond]
        dim = max_bond

    else:
        dim = U.shape[1]

    U = U.reshape(phy_dim, a, dim, order='F')
    V = V.reshape( dim, d, phy_dim, order='F')

    Uprime = np.zeros((phy_dim, dim, a))
    Vprime = np.zeros(( d, dim, phy_dim))

    # making transposed matrices
    for sigma in range(phy_dim):
        Uprime[sigma, :, :] = U[sigma].transpose()
        Vprime[:, :, sigma] = V[:, :, sigma].transpose()

    # explicitly checkying the left and right canonical conditions
    print( 'Left canonical', np.allclose(np.identity(dim), Uprime[0].dot(U[0]) + Uprime[1].dot(U[1]))) 
    print( 'Right canonical',np.allclose(np.identity(dim), V[:, :,0].dot(Vprime[:, :,0]) + V[:, :,1].dot(Vprime[:, :,1]))) 

    # contract U into L
    L = np.einsum( 'ijk,ami,abjl,bkn->mln', L, Uprime, Hmid, U)

    # contract V into R
    R = np.einsum( 'ijk,ima,ablj,nkb->mln', R, Vprime, Hmid, V)

    return L, R, w[0]

# comparing from 10 to 80 from iTensor
def compare():

    ref = [-6.026674183332266, -7.296229810558749, -8.566772233505597, -9.837951447458874, -11.109565585436787, -12.38148999963567, 
    -13.653643543519529, -14.925971109736622, -16.19843395732519, -17.471004053897513, -18.7436606137305, -20.016387897613146, -21.289173768635436, 
    -22.56200871708699, -23.83488518588028, -25.107797093834748, -26.38073949198884, -27.653708311750712, -28.926700177444786, -30.199712264736977, 
    -31.472742192443764, -32.74578793899193, -34.01884777655125, -35.29192022040504, -36.56500398643414, -37.83809795919828, -39.11120116481688, 
    -40.38431274922802, -41.65743196015112, -42.93055813242048, -44.20369067643123, -45.4768290657301, -46.74997283111254, -48.0231215511995, -49.29627484868025, -50.56943238158391]

    L = np.arange(10, 82, 2)
    print(ref/L)

# we run the simplest NN hopping model
if __name__ == '__main__':

    # max bond dim
    max_bond = 20
    # DMRG max iterations
    itrs = 100

    # prepare MPO matrices, with dimensions
    # (2, 2, 4) and (2, 2, 4, 4)
    H = MPOs()

    #L and R states both have dim (2, 2), where 
    #sigmas are at 1 and 2 respectively
    L, R = init_state(H)

    compare()
    for itr in range(itrs):
        L, R, energy = solve(L, R, H, 2, max_bond)
        print(energy)




    