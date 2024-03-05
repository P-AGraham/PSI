import numpy as np 

def factorial(m):
    if m == 0:
        return 1
    return m*factorial(m-1)

def monopole_harmonic(s, m, theta, phi):
    Nm2 = factorial(2 * s + 1)/(4 * np.pi * factorial(s + m) * factorial(s - m))

    c = np.cos(theta/2) ** (s + m)
    s = np.sin(theta/2) ** (s - m)
    
    return np.exp(1j * m * phi) * np.sqrt(Nm2) * c * s