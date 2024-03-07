def creation(m, N):
    c0up = np.kron(np.array([[0, 1], [0, 0]]), np.eye(2))
    c0down = np.kron(np.eye(2), np.array([[0, 1], [0, 0]]))

    cz = np.array([[-1, 0], [0, 1]])

    cup = np.kron(np.kron(np.eye(2**(2*(N-m-1))), c0up), np.eye(2**(2*m)))
    cdown = np.kron(np.kron(np.eye(2**(2*(N-m-1))), c0down), np.eye(2**(2*m)))

    string_up = None
    string_down = None

    for _ in range(0,2*(N-m)):
        string_up = np.kron(string_up, cz)
        string_down = np.kron(string_down, cz)

    string_down = np.kron(string_down, cz)

    string_up = np.kron(string_up, np.eye(2**(2*m)))
    string_down = np.kron(string_down, np.eye(2**(2*m - 1)))
    
    return string_up @ cup, string_down @ cdown

new_bud = cUdag[n2](cD[n3](b))
                new_bdu = cDdag[n2](cU[n3](b))

                new_budud = cUdag[n1](cD[n4](new_bud))
                new_buddu = cDdag[n1](cU[n4](new_bud))

                new_bduud = cUdag[n1](cD[n4](new_bdu))
                new_bdudu = cDdag[n1](cU[n4](new_bdu))

                H_matrix[b, abs(new_budud)] += Vm * new_budud
                H_matrix[b, abs(new_buddu)] += Vm * new_bUD
                H_matrix[b, abs(new_bduud)] += Vm * new_bDU
                H_matrix[b, abs(new_bdudu)] += Vm * new_bDD


#%%
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

print("anihilation")
for i in range(10):
    a1 = c_action(0, 1)[0](i)

    a2 = c_action(0, 1)[0](a1)
    print(bin(i), bin(a1), bin(a2))

print("creation")
for i in range(10):
    a1 = c_action(0, 1)[1](i)

    a2 = c_action(0, 1)[1](a1)
    print(bin(i), bin(a1), bin(a2))

#print(bin(n))
# print(2 * i + s_idx)
# print(bin(n)[:-(2 * i + s_idx)])
# ask about local and ultra local, ask about the coupling depending on the number of electrons


def factorial(m):
    if m == 0:
        return 1
    return m*factorial(m-1)

def wigner_3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0:
        return 0
    
    elif np.abs(j1-j2) > j3:
        return 0
    
    elif j3 > j1 + j2: 
        return 0
    
    f1 = factorial(j1 + j2 - j3) * factorial(j1 - j2 + j3) * factorial(-j1 + j2 + j3) 
    f2 = factorial(j1 + j2 + j3 + 1)
    f31 = factorial(j1 - m1) * factorial(j1 + m1)
    f32 = factorial(j2 - m2) * factorial(j2 + m2)
    f33 = factorial(j3 - m3) * factorial(j3 + m3)
    f3 = f31 * f32 * f33

    S  = 0

    for k in range(0, N + 1):
        f4 = factorial(j1 + j2 - j3 -k) 
        f5 = factorial(j1 - m1 - k) * factorial(j2 + m2 - k) 
        f6 = factorial(j3 - j2 + m1 + k) * factorial(j3 - j1 -  m2 + k)
        S += (-1) ** k /(factorial(k) * np.sqrt(f4 * f5 * f6))

    return (-1) ** (j1 - j2 - m3) * np.sqrt(f1/f2) * np.sqrt(f3) * S


 new_i = 0
    for i in range(dim):
        
        new_j = 0
        if bin(i).count('1') == N:
            for j in range(dim):
                if bin(j).count('1') == N:
                    H_hfill[new_i, new_j] += H_matrix[i, j]
                    
                    new_j += 1
            new_i += 1

def test_sector(n, N):
    n_bin = bin(n)[2:]
    test = True
    for i in range(n//2):
        test_i = n_bin[2*i:2*i+1].count('1') <= 1
        test = test and test_i
    
    return test and n_bin.count('1')


bin_i = bin(i)[2:][::-1] 
            bin_i = bin_i + '0' * (2 * N - len(bin_i))
            states_label.append(bin_i)


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
        Z_srep[b, b] = s - b

    for b in range(int(2*s)):
        X_srep[b, b + 1] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/2
        X_srep[b + 1, b] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/2
        Y_srep[b, b + 1] = -np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/(2j)
        Y_srep[b + 1, b] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/(2j)


    Pair = zip([X_srep, Y_srep, Z_srep], [Lx_matrix, Ly_matrix, Lz_matrix])

    for axis, matrix in Pair:
        for m1 in m:
            for m2 in m:
                n1, n2 = int(m1+s), int(m2+s)
                for b in range(dim):

                    new_bDD, sign_bDD = cDdag[n2](*cD[n1](b, 1))
                    matrix[b, new_bDD] += axis[n1, n2] * sign_bDD

                    new_bUU, sign_bUU = cUdag[n2](*cU[n1](b, 1))
                    matrix[b, new_bUU] += axis[n1, n2] * sign_bUU

    Lx_hfill, _ = hfill_sector(Lx_matrix, N)
    Ly_hfill, _ = hfill_sector(Ly_matrix, N)
    Lz_hfill, _ = hfill_sector(Lz_matrix, N)

    return Lx_hfill, Ly_hfill, Lz_hfill


count = 0
for i in range(len(L2)):
    #print(L2[i, i])
    if (L2[i, i]-6)**2 < 1:
        A = L2[i, i]
        J = np.real(j(A))
        print(L2[i, i], J, "\t", np.real(Lz_hfill[i, i]))
        count += 1

print(count)


# %%
N = 5
s = (N - 1)/2
dim = 2 ** (2 * N) 

Lx_matrix = np.zeros((dim, dim), dtype=complex)
Ly_matrix = np.zeros((dim, dim), dtype=complex)
Lz_matrix = np.zeros((dim, dim), dtype=complex)

m = np.arange(-s, s+1, 1)

X_srep = np.zeros((int(2*s+1), int(2*s+1)), dtype=complex)
Y_srep = np.zeros((int(2*s+1), int(2*s+1)), dtype=complex)
Z_srep = np.zeros((int(2*s+1), int(2*s+1)), dtype=complex)

for b in range(int(2*s+1)):
    Z_srep[b, b] = s - b

for b in range(int(2*s)):
    X_srep[b, b + 1] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/2
    X_srep[b + 1, b] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/2
    Y_srep[b, b + 1] = -np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/(2j)
    Y_srep[b + 1, b] =  np.sqrt(2 * (s + 1) * (b + 1) - (b + 2) * (b + 1))/(2j)

print(X_srep)
print(Y_srep)
print(Z_srep)
print(X_srep @ X_srep + Y_srep @ Y_srep + Z_srep @ Z_srep)


Any spherically symmetric quantity can be expanded in spherical harmonics or monopole harmonics. For a given $l$, we have non-zero monopole harmonics for $s \in [0, l]$ and their conjugate pair (with negative $l$ in the spherical harmonics used to construct them. I suspect there is a simpler formulation of this adapted to eq. 3 of the paper) so $2s + 1$ times more ... It appears that the combined spherical harmonics in $l\in [0, 2s]$ span the same space as the monopole harmonics with $l = s$ (there are $2s$ of them)  