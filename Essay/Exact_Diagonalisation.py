#%%

import numpy as np 
from sympy.physics.wigner import wigner_3j
from scipy.sparse.linalg import eigs
from scipy.special import comb

s = 1 # Half integer, 2s is the max landau level
N = 4 # number of electrons
s = (N - 1)/2


V0 = 4.75
V1 = 1


def V(m1, m2, m3, m4, V0, V1 = 1):

    W10 = wigner_3j(s, s, 2 * s - 0, m1, m2 - m1, - m1 - m2)
    W20 = wigner_3j(s, s, 2 * s - 0, m4, m3, - m3 - m4)
    W11 = wigner_3j(s, s, 2 * s - 1, m1, m2 - m1, - m1 - m2)
    W21 = wigner_3j(s, s, 2 * s - 1, m4, m3, - m3 - m4)

    f0 = 4 * s - 2 * 0 + 1
    f1 = 4 * s - 2 * 1 + 1

    return f0 * V0 * W10 * W20 + f1 *  V1 * W11 * W21


def c_action(i, s):
    # s = 0, 1

    def anihilate_site_spin(n):

        if n >= 0:
            n_bin = bin(n)[2:][::-1]
            sign = 1
        else:
            n_bin = bin(-n)[2:][::-1]
            sign = -1

        if 2 * i + s < len(n_bin):
            filled = n_bin[2 * i + s]
        else: 
            filled = 0
        
        if not int(filled):
            new_n = 0 
        
        else:
            anti_c = n_bin[:2 * i + s].count('1')
            new_n = n - sign * 2 ** (2 * i + s)
            new_n = new_n * (-1) ** anti_c
        
        return new_n 
    
    def create_site_spin(n):

        if n >= 0:
            n_bin = bin(n)[2:][::-1]
            sign = 1
        else:
            n_bin = bin(-n)[2:][::-1]
            sign = -1

        if 2 * i + s < len(n_bin):
            filled = n_bin[2 * i + s]
        else: 
            filled = 0
        
        if int(filled):
            new_n = 0 
        
        else:
            anti_c = n_bin[:2 * i + s].count('1')
            new_n = n + sign * 2 ** (2 * i + s)
            new_n = new_n * (-1) ** anti_c
        
        return new_n
    
    return anihilate_site_spin,  create_site_spin



def H(h, N, V0, V1 = 1):
    s = (N - 1)/2
    dim = 2 ** (2 * N) 
    hfill_dim = int(comb(2 * N, N))

    H_matrix = np.zeros((dim, dim))
    H_hfill = np.zeros((hfill_dim, hfill_dim)) 

    m = np.arange(-s, s+1, 1)
    m_arr = np.array(np.meshgrid(*[m, m, m, m])).T.reshape(-1, 4)


    cU = []
    cD = []
    cUdag = []
    cDdag = []

    for n in range(0, N):
        
        c, d = c_action(n, 0)
        cD.append(c) 
        cDdag.append(d)

        a, b = c_action(n, 1)
        cU.append(a)
        cUdag.append(b)


    for m1, m2, m3, m4 in m_arr:
        n1, n2, n3, n4 = int(m1+s), int(m2+s), int(m3+s), int(m4+s)

        if m2 == m3 == m4 == 0:
            for b in range(dim):
                new_bUD = cUdag[n1](cD[n1](b))
                H_matrix[b, abs(new_bUD)] += -h * new_bUD

                new_bDU = cDdag[n1](cU[n1](b))
                H_matrix[b, abs(new_bDU)] += -h * new_bDU

        if m1 + m2 == m3 + m4: 

            Vm = V(m1, m2, m3, m4, V0, V1 = V1)/2

            for b in range(dim):
                
                # H00
                
                new_bU = cUdag[n2](cU[n3](b))
                new_bD = cDdag[n2](cD[n3](b))

                new_bUU = cUdag[n1](cU[n4](new_bU))
                new_bUD = cDdag[n1](cD[n4](new_bU))

                new_bDU = cUdag[n1](cU[n4](new_bD))
                new_bDD = cDdag[n1](cD[n4](new_bD))

                #H_matrix[b, abs(new_bUU)] += Vm * new_bUU
                #H_matrix[b, abs(new_bUD)] += Vm * new_bUD
                #H_matrix[b, abs(new_bDU)] += Vm * new_bDU
                #H_matrix[b, abs(new_bDD)] += Vm * new_bDD
#
                ## Hzz (dont forget the sign)
#
                #H_matrix[b, abs(new_bUU)] += (1) * (1) * Vm * new_bUU
                #H_matrix[b, abs(new_bUD)] += (1) * (-1) * Vm * new_bUD
                #H_matrix[b, abs(new_bDU)] += (-1) * (1) * Vm * new_bDU
                #H_matrix[b, abs(new_bDD)] += (-1) * (-1) * Vm * new_bDD

    print(H_matrix)
    new_i = 0
    for i in range(dim):
        
        new_j = 0
        if bin(i).count('1') == N:
            for j in range(dim):
                if bin(j).count('1') == N:
                    H_hfill[new_i, new_j] += H_matrix[i, j]
                    
                    new_j += 1
            new_i += 1

    return H_hfill


H_solve = H(3, 1, 4.75, 1)
print(H_solve)
#%%
print(np.round(eigs(H_solve, which='SM')[1].T, 5))

# %%
def factorial(m):
    if m == 0:
        return 1
    return m*factorial(m-1)

factorial(3)
# %%
