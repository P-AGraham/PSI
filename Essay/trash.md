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