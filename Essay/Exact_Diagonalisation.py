import numpy as np 

s = 1 # Half integer, 2s is the max landau level
N = 4 # number of electrons
s = (N - 1)/2


V0 = 4.75
V1 = 1


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



def V(m1, m2, m3, m4, V0, V1 = 1):

    W10 = wigner_3j(s, s, 2 * s - 0, m1, m2 - m1, - m1 - m2)
    W20 = wigner_3j(s, s, 2 * s - 0, m4, m3, - m3 - m4)
    W11 = wigner_3j(s, s, 2 * s - 1, m1, m2 - m1, - m1 - m2)
    W21 = wigner_3j(s, s, 2 * s - 1, m4, m3, - m3 - m4)

    return V0 * W10 * W20 + V1 * W11 * W21


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




def H_transverse(h, N, V0, V1 = 1):
    s = (N - 1)/2
    dim = 2 ** (2 * N)
    H_transverse = np.zeros((dim, dim))
    H_densitydensity = np.zeros((dim, dim))
    H_spinspin = np.zeros((dim, dim))


    m = np.arange(-s, s+1, 1)
    m_arr = np.array(np.meshgrid(*[m, m, m, m])).T.reshape(-1, 4)

    for m1, m2, m3, m4 in m_arr:

        c1up, c1down = creation(m1, N)
        c2up, c2down = creation(m2, N)
        c3up, c3down = creation(m3, N)
        c4up, c4down = creation(m4, N)


        if m2 == m3 == m4 == 0:
            H_transverse += -h * (c1up.T * c1down + c1down.T * c1up)

        if m1 + m2 == m3 + m4: 
            Vm = V(m1, m2, m3, m4, V0, V1 = V1)
            c14_up = c1up.T @ c4up
            c14_down = c1down.T @ c4down
            c23_up = c2up.T @ c3up
            c23_down = c2down.T @ c3down



            H_densitydensity += Vm *  @ (c1up.T @ c4up + c1down.T @ c4down)

            H_spinspin += Vm * (c1.T @ sigma_z @ c4) @ (c2.T @ sigma_z @ c3)

    return H




#%%
import numpy as np

print(1)

#c1 @ np.array([0, 0, 1, 0])
# %%


# ask about local and ultra local, ask about the coupling depending on the number of electrons
# %%
