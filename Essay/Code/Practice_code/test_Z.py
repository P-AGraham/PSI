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
            new_bin_b += bin_b[i+1]
            new_bin_b += bin_b[i]
            sign += (int(bin_b[i]) * int(bin_b[i+1]))
        #print(sign)

        new_b = int("0b"+new_bin_b, 2)
        Z2_matrix[b, new_b] += (-1) ** sign

    Z2_hfill, _ = hfill_sector(Z2_matrix, N)

    return Z2_hfill


def Z2(N):
    s = (N - 1)/2
    dim = 2 ** (2 * N) 

    Z2_matrix = np.zeros((dim, dim), dtype=complex)

    m = np.arange(-s, s+1, 1)

    cU, cD, cUdag, cDdag = c_wrap_list(N)

    for b in range(dim):
        b1 = b
        for m1 in m:
            n1 = int(s+m1)

            new_bDU, sign_bDU = cUdag[n1](*cD[n1](b1, 1))
            new_bUD, sign_bUD = cDdag[n1](*cU[n1](b1, 1))

            Z2_matrix[b, new_bDU] += sign_bDU
            Z2_matrix[b, new_bUD] += sign_bUD

    Z2_hfill, _ = hfill_sector(Z2_matrix, N)

    return Z2_hfill