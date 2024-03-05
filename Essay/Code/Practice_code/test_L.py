#%%
import numpy as np 
from sympy.physics.wigner import wigner_3j
from scipy.sparse.linalg import eigs
from scipy.special import comb


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


def L(N):
    s = (N - 1)/2
    dim = 2 ** (2 * N) 

    Lx_matrix = np.zeros((dim, dim), dtype=complex)
    Ly_matrix = np.zeros((dim, dim), dtype=complex)
    Lz_matrix = np.zeros((dim, dim), dtype=complex)

    Lplus_matrix = np.zeros((dim, dim), dtype=complex)
    Lz2_matrix = np.zeros((dim, dim), dtype=complex)

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

    for m1 in m[:-1]:
        n1 = int(s + m1)

        for b in range(dim):
            new_bDD, sign_bDD = cDdag[n1 + 1](*cD[n1](b, 1))
            Lplus_matrix[b, new_bDD] += np.sqrt((s-m1)*(s+m1+1)) * sign_bDD

            new_bUU, sign_bUU = cUdag[n1 + 1](*cU[n1](b, 1))
            Lplus_matrix[b, new_bUU] += np.sqrt((s-m1)*(s+m1+1)) * sign_bUU

    for m1 in m:
        n1 = int(s + m1)

        for b in range(dim):
            new_bDD, sign_bDD = cDdag[n1](*cD[n1](b, 1))
            Lz2_matrix[b, new_bDD] += m1 * sign_bDD

            new_bUU, sign_bUU = cUdag[n1](*cU[n1](b, 1))
            Lz2_matrix[b, new_bUU] += m1 * sign_bUU
            
    Lx_hfill, _ = Lx_matrix, 0#hfill_sector(Lx_matrix, N)
    Ly_hfill, _ = Ly_matrix, 0#hfill_sector(Ly_matrix, N)
    Lz_hfill, _ = Lz_matrix, 0#hfill_sector(Lz_matrix, N)

    return Lx_hfill, Ly_hfill, Lz_hfill, Lplus_matrix, Lz2_matrix

N = 2
s = (N - 1)/2
dim = 2 ** (2 * N) 

Lx_hfill, Ly_hfill, Lz_hfill, Lplus_matrix, Lz2_matrix = L(N)
L2 = Lx_hfill @ Lx_hfill + Ly_hfill @ Ly_hfill + Lz_hfill @ Lz_hfill
L2_2 =  Lz2_matrix @ (Lz2_matrix + np.eye(dim)) + 2 * Lplus_matrix @ Lplus_matrix.T


L2, _ = hfill_sector(L2, N)
L2_2, _ = hfill_sector(L2_2, N)

Lx_hfill, _ = hfill_sector(Lx_hfill, N)
Ly_hfill, _ = hfill_sector(Ly_hfill, N)
Lz_hfill, state_lbl = hfill_sector(Lz_hfill, N)

print(np.all(np.round(Lx_hfill @ Ly_hfill - Ly_hfill @ Lx_hfill, 3) == 1j * np.round(Lz_hfill, 3)))
print(np.all(np.round(Ly_hfill @ Lz_hfill - Lz_hfill @ Ly_hfill, 3) == 1j * np.round(Lx_hfill, 3)))
print(np.all(np.round(Lz_hfill @ Lx_hfill - Lx_hfill @ Lz_hfill, 3) == 1j * np.round(Ly_hfill, 3)))

print(np.all(np.round(Lz_hfill @ L2 - L2 @ Lz_hfill, 2) == 0))
#print(np.all(np.round(Lz_hfill @ L2_2 - L2_2 @ Lz_hfill, 2) == 0))

print(np.all(np.round(Ly_hfill @ L2 - L2 @ Ly_hfill, 2) == 0))
#print(np.all(np.round(Ly_hfill @ L2_2 - L2_2 @ Ly_hfill, 2) == 0))

print(np.all(np.round(Lx_hfill @ L2 - L2 @ Lx_hfill, 2) == 0))
#print(np.all(np.round(Lx_hfill @ L2_2 - L2_2 @ Lx_hfill, 2) == 0))

def j(j2j):
    return (-1 + np.sqrt(1 + 4 * j2j))/2

for i in range(len(L2)):
    print(L2[i,i], Lz_hfill[i, i])
print([bin(n)[2:][::-1] for n in range(dim)])
# %%
np.all(np.round(L2 - np.diag(np.diag(L2)), 3) == 0)
# %%
eig_val, eig_vec = eigs(L2, 70, which='SR')
eig_val = np.real(np.round(eig_val, 6))
eig_val_sort = np.sort(eig_val)


idx = eig_val.argsort()
eig_vec_sort = eig_vec[idx]

eig_val_sort_Lz = np.diag(np.dot(np.dot(eig_vec_sort.T, Lz_hfill), (eig_vec_sort)))

print(eig_val_sort, eig_val_sort_Lz)
# %%
