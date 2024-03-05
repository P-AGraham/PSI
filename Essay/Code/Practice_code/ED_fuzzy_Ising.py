#%%
import numpy as np 
from sympy.physics.wigner import wigner_3j
from scipy.sparse.linalg import eigsh
from scipy.special import comb
import pandas as pd 
from IPython.display import display


def V(m1, m2, m3, m4, s, V0, V1 = 1):

    W10 = wigner_3j(s, s, 2 * s - 0, m1, m2, - m1 - m2)
    W20 = wigner_3j(s, s, 2 * s - 0, m4, m3, - m3 - m4)
    W11 = wigner_3j(s, s, 2 * s - 1, m1, m2, - m1 - m2)
    W21 = wigner_3j(s, s, 2 * s - 1, m4, m3, - m3 - m4)

    f0 = 4 * s - 2 * 0 + 1
    f1 = 4 * s - 2 * 1 + 1

    return f0 * V0 * W10 * W20 + f1 *  V1 * W11 * W21


def c_wrap(i, spin): 
    # s : spin is 0 for down and 1 for up 
    # i : site where c acts from 0 to (2*s+1 - 1)
    # (i, s) gives the location of the target 
    # n : state on which c acts 
    # sign : coefficient of the state

    def c(n, sign):
    
        if sign == 0: return n, sign 
         # if the state is the null state, there is nothing to do

        n_bin = bin(n)[2:][::-1]
        # binary expression of n, reversed to keep operator action convention

        if 2 * i + spin < len(n_bin): filled = n_bin[2 * i + spin]
        else: filled = 0
        # if the target is out of binary expansion, return null
    
        if not int(filled): return n, 0
    
        else:
            # count the number of "1" before the target
            anti_c = n_bin[:2 * i + spin].count('1')
            # Removes an electron at the filled target
            new_n = n - 2 ** (2 * i + spin)
            return new_n, sign * (-1) ** anti_c
    
    def cdag(n, sign):
    
        if sign == 0: return n, sign 

        n_bin = bin(n)[2:][::-1]

        if 2 * i + spin < len(n_bin): filled = n_bin[2 * i + spin]
        else: filled = 0
    
        if int(filled): return n, 0
        else:
            anti_c = n_bin[:2 * i + spin].count('1')
            new_n = n + 2 ** (2 * i + spin)
            return new_n, sign * (-1) ** anti_c
        
    return c, cdag

def c_wrap_list(N):
    cU, cD, cUdag, cDdag = [], [], [], []

    for i in range(0, N):
    
        a, b = c_wrap(i, 1)
        cU.append(a)
        cUdag.append(b)

        c, d = c_wrap(i, 0)
        cD.append(c) 
        cDdag.append(d)

    return cU, cD, cUdag, cDdag


def hfill_sector(matrix, N):
    dim = 2 ** (2 * N) 
    hfill_dim = int(comb(2 * N, N))

    matrix_hfill = np.zeros((hfill_dim, hfill_dim), dtype=complex)

    new_i = 0
    states_label = []
    for i in range(dim):
        new_j = 0
        
        if bin(i).count('1') == N:
            bin_i = bin(i)[2:][::-1] 
            bin_i = bin_i + '0' * (2 * N - len(bin_i))
            states_label.append(bin_i)
            
            for j in range(dim):
                if bin(j).count('1') == N:
                    matrix_hfill[new_i, new_j] += matrix[i, j]
                    
                    new_j += 1
            new_i += 1

    return matrix_hfill, states_label

import numpy as np

def Z2(N):
    dim = 2 ** (2 * N)

    Z2_matrix = np.zeros((dim, dim), dtype=complex)

    for b in range(dim):
        bin_b = bin(b)[2:]
        #print(bin_b)
        bin_b = '0' * (2 * N - len(bin_b)) + bin_b
        new_bin_b = ""

        sign = 0
        for i in range(0, 2 * N - 1, 2):
            # Construct the string with spin up and down on site i exchanged
            new_bin_b += bin_b[i+1]
            new_bin_b += bin_b[i]
            # Get a minus if the swaped states are both occupied
            sign += (int(bin_b[i]) * int(bin_b[i+1])) 

        new_b = int("0b"+new_bin_b, 2)
        Z2_matrix[b, new_b] += (-1) ** sign

    Z2_hfill, _ = hfill_sector(Z2_matrix, N)

    return Z2_hfill

def L(N):
    s = (N - 1)/2
    dim = 2 ** (2 * N) 

    Lx_matrix = np.zeros((dim, dim), dtype=complex)
    Ly_matrix = np.zeros((dim, dim), dtype=complex)
    Lz_matrix = np.zeros((dim, dim), dtype=complex)

    m = np.arange(-s, s+1, 1)

    cU, cD, cUdag, cDdag = c_wrap_list(N)

    X_srep = np.zeros((int(2*s+1), int(2*s+1)), dtype=complex)
    Y_srep = np.zeros((int(2*s+1), int(2*s+1)), dtype=complex)
    Z_srep = np.zeros((int(2*s+1), int(2*s+1)), dtype=complex)

    for b in range(int(2*s+1)):
        Z_srep[b, b] = b - s

    for b in range(int(2*s)):
        X_srep[b, b + 1] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/2
        X_srep[b + 1, b] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/2
        Y_srep[b, b + 1] = -np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/(2j)
        Y_srep[b + 1, b] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/(2j)


    Pair = zip([X_srep, Y_srep, Z_srep], [Lx_matrix, Ly_matrix, Lz_matrix])

    for axis, matrix in Pair:
        for m1 in m:
            for m2 in m:
                n1, n2 = int(s+m1), int(s+m2)

                for b in range(dim):

                    new_bDD, sign_bDD = cDdag[n2](*cD[n1](b, 1))
                    matrix[b, new_bDD] += axis[n1, n2] * sign_bDD

                    new_bUU, sign_bUU = cUdag[n2](*cU[n1](b, 1))
                    matrix[b, new_bUU] += axis[n1, n2] * sign_bUU

    Lx_hfill, _ = hfill_sector(Lx_matrix, N)
    Ly_hfill, _ = hfill_sector(Ly_matrix, N)
    Lz_hfill, _ = hfill_sector(Lz_matrix, N)

    return Lx_hfill, Ly_hfill, Lz_hfill


def H(N, h, V0, V1 = 1):
    s = (N - 1)/2
    dim = 2 ** (2 * N) 
    hfill_dim = int(comb(2 * N, N))

    H_matrix = np.zeros((dim, dim))
    H_hfill = np.zeros((hfill_dim, hfill_dim)) 

    m = np.arange(-s, s+1, 1)
    m_arr = np.array(np.meshgrid(*[m, m, m, m])).T.reshape(-1, 4)

    cU, cD, cUdag, cDdag = c_wrap_list(N)

    for m1 in m:
        n1 = int(m1+s)
        for b in range(dim):
            new_bUD, sign_bUD = cUdag[n1](*cD[n1](b, 1))
            H_matrix[b, new_bUD] += -h * sign_bUD

            new_bDU, sign_bDU = cDdag[n1](*cU[n1](b, 1))
            H_matrix[b, new_bDU] += -h * sign_bDU

    for m1, m2, m3, m4 in m_arr:
        if m1 + m2 == m3 + m4: 

            n1, n2, n3, n4 = int(m1+s), int(m2+s), int(m3+s), int(m4+s)
            Vm = V(m1, m2, m3, m4, s, V0, V1 = V1)/2

            for b in range(dim):
                
                new_bU = cUdag[n2](*cU[n3](b, 1))
                new_bD = cDdag[n2](*cD[n3](b, 1))

                new_bUU, sign_bUU = cUdag[n1](*cU[n4](*new_bU))
                new_bUD, sign_bUD = cDdag[n1](*cD[n4](*new_bU))

                new_bDU, sign_bDU = cUdag[n1](*cU[n4](*new_bD))
                new_bDD, sign_bDD = cDdag[n1](*cD[n4](*new_bD))

                # H00
                H_matrix[b, new_bUU] += Vm * sign_bUU
                H_matrix[b, new_bUD] += Vm * sign_bUD
                H_matrix[b, new_bDU] += Vm * sign_bDU
                H_matrix[b, new_bDD] += Vm * sign_bDD

                ## Hzz (dont forget the sign)
                H_matrix[b, new_bUU] += -( 1) * ( 1) * Vm * sign_bUU
                H_matrix[b, new_bUD] += -( 1) * (-1) * Vm * sign_bUD
                H_matrix[b, new_bDU] += -(-1) * ( 1) * Vm * sign_bDU
                H_matrix[b, new_bDD] += -(-1) * (-1) * Vm * sign_bDD

    H_hfill, states_label = hfill_sector(H_matrix, N)

    return H_hfill, states_label


def solve(N, h, V0, V1 = 1, k = 70):
    print(locals())
    
    # Creating the Hamiltonian matrix in the half filling sector
    H_hfill, _ = H(N, h, V0, V1 = V1)

    v0 = np.random.rand(int(comb(2 * N, N)))

    # Sparse Exact Diagonalization
    eig_val, eig_vec = eigsh(H_hfill, k = k, which='SR', v0=v0)

    
    
    # Eigen value sorting index
    eig_val = np.real(eig_val)
    idx = eig_val.argsort()

    # Sorting eigen vectors and eigen values with index
    eig_val_sort = eig_val[idx]
    eig_vec_sort = eig_vec.T[idx]

    # Casimir operator for SO(3) symmetry 
    Lx_hfill, Ly_hfill, Lz_hfill = L(N)
    L2 = Lx_hfill @ Lx_hfill + Ly_hfill @ Ly_hfill + Lz_hfill @ Lz_hfill

    ell2_ell = [(np.conj(vec).T @ L2 @ vec)/(np.conj(vec).T @ vec) for vec in eig_vec_sort]

    ell = np.round((-1 + np.sqrt(1 + 4 * np.array(ell2_ell)))/2, 2)
    ell = np.real(ell)

    # Ising Z2 symmetry

    Z2_hfill = Z2(N)
    Z2_val = np.zeros(k)

    i = 0
    while i < k: 
        
        stop_i = int(2*ell[i]+1)
        print(i, i +  stop_i)

        cluster = eig_vec_sort[i:i+stop_i]
        cluster = np.array([c/np.sqrt(np.dot(np.conj(c), c)) for c in cluster])
    
        Z2_cluster = np.conj(cluster) @ Z2_hfill @ cluster.T
        Z2_val_cluster, _ = np.linalg.eigh(Z2_cluster)
        Z2_val[i:i+stop_i] = Z2_val_cluster
        
        i += stop_i

    Z2_val = np.round(Z2_val, 4)

    print(np.all(Z2_hfill.T == Z2_hfill))
    #print(np.round(eigsh(Z2_hfill, k = k, which='SR', v0=v0)[0]))
    #print(np.all(np.round(Z2_hfill @ L2 - L2 @ Z2_hfill, 3)==0))
    #print(np.all(np.round(Z2_hfill @ H_hfill - H_hfill @ Z2_hfill, 3)==0))
    #print(np.all(np.round(Z2_hfill @ Lz_hfill - Lz_hfill @ Z2_hfill, 3)==0))
    # Particle-hole symmetry
    #P = [0]*k

    # Scaling dimension (Rescale the spectrum shifts so that the SET energy gap is 3 (3D CFT))
        
    Delta = np.zeros(k)

    #SET_loc = np.where((ell == 2) & (Z2_val > 0))[0][0]
    #E_shift = [eig_val_sort[i+1] - eig_val_sort[i] for i in range(k-1)]
    #SET_E_shift = E_shift[SET_loc]

    #Delta[:-1] = E_shift/SET_E_shift*3
    

    table = pd.DataFrame({
        "Energy"           : eig_val_sort, 
        "Orbital Momentum" : ell, 
        "Ising Z2"         : Z2_val, 
        #"Particle-Hole"    : P,
        "Scaling Dimension": Delta,
        })
    
    return table

if __name__ == "__main__":
    table = solve(4, 3.16, 4.75, k = 50)
    display(table)
    
# %%
