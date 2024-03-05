#%%
import numpy as np
from scipy.special import jv
import matplotlib.pyplot as plt

def I_m(k, a, T, omega_1, omega_2, omega_3, n_cutoff = 10):

    tau_i = np.array([-T, -T/2, T/2, T])
    omega_i = np.array([omega_1, omega_2, omega_3])
    gamma_i = 1/np.sqrt(1-omega_i**2 * a**2)

    delta_i = [np.sum(tau_i[1:i+1]*omega_i[0:i]) for i in range(1, len(omega_i))]
    delta_i = np.array([0]+delta_i)

    def Ifunc(Omega):
        I = 0
        Omega_i = gamma_i * k + Omega
        for n in range(-n_cutoff, n_cutoff+1):

            Jn = jv(n, -k * a)

            E1 = np.exp(1j * n * delta_i)

            E2_m = np.exp(1j * (Omega_i + n * omega_i) * tau_i[0:-1])
            E2_p = np.exp(1j * (Omega_i + n * omega_i) * tau_i[1:])

            denom = 1j * (Omega_i + n * omega_i)
            print(denom)
            print(Omega_i)

            I += Jn * np.sum(E1 * (E2_p - E2_m)/denom) * (1j) ** n
        return I 
    
    return Ifunc


a = 0.5
If = I_m(1, a, 10, 1.9, 0, 1.9)
Omega_lin = np.linspace(1, 5, 200)

plt.plot(Omega_lin, [np.real(If(O)) for O in Omega_lin], "r")
plt.plot(Omega_lin, [np.imag(If(O)) for O in Omega_lin], "b")

#plt.plot(Omega_lin, [np.sqrt(np.imag(If(O))**2 + np.real(If(O))**2) for O in Omega_lin])

# %%
