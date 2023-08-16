import matplotlib.pyplot as plt
import matplotlib as mpl

def gen_op(L, onesiteop, twositeop):

    if len(onesiteop) > 1:
        print('one site op # error')
        exit()

    state = 2

    #first pass, get the total state space
    for _, _, _, int_range in twositeop:
        state += int_range

    A = [[['0'] for _ in range(state)] for _ in range(state)]
    start = [[['0'] for _ in range(state)]]
    end = [[['0']] for _ in range(state)]

    # second pass, fill the matrices

    # set I
    A[0][0] = ['I']
    end[0][0] = ['I' + '_{{{}}}'.format(L)]
    start[0][-1] = ['I']
    A[-1][-1] = ['I']

    cnt = 1

    for op1, op2, int_st, int_range in twositeop:


        if len(int_st) != int_range:
            raise Exception('not enough int parameters for int_range')

        A[cnt][0] = [op1]
        end[cnt][0] = [op1 + '_{{{}}}'.format(L)]
        A[-1][cnt] = start[0][cnt] = [int_st[0] + op2]
        
        cnt += 1

        for i in range(int_range-1):
            A[cnt][0] = end[cnt][0] = ['0']
            A[-1][cnt] = start[0][cnt] = [ int_st[i + 1] + op2]

            A[cnt ][cnt - 1 ] = ['I']

            cnt += 1


    # set one site operator
    A[-1][0] =  [ onesiteop[0]]
    start[0][0] = [ onesiteop[0]]
    end[-1][0] = [ onesiteop[0] + '_{{{}}}'.format(L)]

    #print(start, A, end)
    return start, A, end

def mul(A, B, L, idflag):

    #A: n * k
    #B: k * m

    n = len(A)
    k = len(B)
    m = len(B[0])

    if k != len(A[0]):

        print('dimensions mismatch')
        exit()

    res = [[[] for _ in range(m)] for _ in range(n)] 

    for i in range(n):
        for j in range(m):

            for p in range(k):
                
                for s in A[i][p]:
                    
                    if s[0] == '0':
                        res[i][j].append('0')

                    else:
                        for t in B[p][j]:
                            if t[0] == '0':
                                res[i][j].append('0')

                            elif not idflag and t[0] == 'I':
                                res[i][j].append(s + '_{{{}}}'.format(L))

                            else:

                                if not idflag and s[0] == 'I':
                                    res[i][j].append(t)

                                else:
                                    res[i][j].append(s  + '_{{{}}}'.format(L) +  t)

    return res
    

def tostring(arr):
    arr = arr[0][0]
    res = ''
    cnt = 1

    for it in arr:

        if len(res) // 193 > cnt:
            cnt += 1
            res += ' $\n$ '
        
        if it != '0':

            if not res:
                res = 'H = ' + it

            else:
                res += '+' + it

    return res


def init(L):

    #add onesite, can only be one!
    onesiteop = ['h']
    #onesiteop = ['0']

    # add two-site operator, following the format: op1, op2, array of interaction strength, range of the interaction
    twositeop = [ ['c', 'c^{{\dagger}}', ['t'], 1], \
                    ['c^{{\dagger}}', 'c', ['t'], 1], \
                        ['\hat{{n}}', '\hat{{n}}', [r'\frac{{1}}{{r_{{{}}}}}'.format(i) for i in range(1, L )], L - 1]]
    start, A, end = gen_op(L, onesiteop, twositeop)


    return start, A, end

def construct(start, A, end, L, idflag):

    if L < 2:
        print('cannot construct 2-site')
        exit()

    cur = mul(A, end, L, idflag)
    L -= 1

    while L > 1:
        cur = mul(A, cur, L, idflag)
        L -= 1

    cur = mul(start, cur, L, idflag)
    return cur



def plotall(st, start, A, end):

    print(st)

    mpl.rcParams['text.latex.preamble'] = r"\usepackage{{amsmath}}"
    mpl.rcParams['text.usetex'] = True

    fig, ax = plt.subplots()
    ax.text(0, 0.5, "${}$".format(st))
    ax.axis('off')

    # convert python array to latex string
    def convert(mat):

        res = r'$\begin{pmatrix}'
        
        for row in mat:
            for i, col in enumerate(row):

                for it in col:

                    res += ' ' + it 

                    if i < len(row) - 1:
                        res += ' & '

            if len(mat) > 1:
                res += r' \\ '
        
        res += r'\end{pmatrix}$'
        return res

    print("start", convert(start))
    print("A", convert(A))
    print("end", convert(end))

    ax.text(0, 0.3, 'Start = {}, A = {}, end = {}'.format(convert(start), convert(A), convert(end)))
    #ax.text(0, 0.3, r'$\begin{pmatrix} ' + r'1 & 2 & 3 \end{pmatrix}$')
    plt.show()

if __name__ == '__main__':

    L = 7

    start, A, end = init(L)
    idflag = True
    res = construct(start, A, end, L - 1, idflag)

    st = tostring(res)
    plotall(st, start, A, end)

    
    #print(mul(start, mul(A,end)))
