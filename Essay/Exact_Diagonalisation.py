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


def Pseudo_V(l):
    pass 
    #Vshort_range

def V(m1, m2, m3, m4):
    for l in range():
        W1 = wigner_3j(s, s, 2 * s - l, m1, m2 - m1, - m1 - m2)
        W2 = wigner_3j(s, s, 2 * s - l, m4, m3, - m3 - m4)
    return



def creation(m, N):
    return 

def H_transverse():
    return 

def H_densitydensity():
    return 

def H_densitydensity():
    return 

#%%
import numpy as np


def creation(m, N):
    c0 = np.array([[0, 1], [0, 0]])
    cz = np.array([[-1, 0], [0, 1]])

    c = np.kron(np.kron(np.eye(2**(N-m-1)), c0), np.eye(2**m))

    string = 1

    for i in range(0,N-m):
        string = np.kron(string, cz)
        #print(string) 

    string = np.kron(string, np.eye(2**m))
    
    return string @ c


c0 = creation(0, 10).T
c1 = creation(1, 10).T

print(c0 @ c1 + c1 @ c0)

#c1 @ np.array([0, 0, 1, 0])
# %%
# ask about local and ultra local, ask about the coupling depending on the number of electrons