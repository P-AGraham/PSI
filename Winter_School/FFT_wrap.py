#%%
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import quad 


def peicewise(Xi, Ti):
    # Xi and Ti are the lab coordinates at each turning point

    def traj(tau):

        x = np.nan
        t = np.nan

        # acumulated propertime
        taui = 0
        
        for i, T in enumerate(Ti[:-1]):

            T_start = T 
            T_end = Ti[i + 1]

            X_start = Xi[i]
            X_end = Xi[i + 1]

            V = (X_end - X_start)/(T_end-T_start)

            G = 1/np.sqrt(1 - np.dot(V, V))
            

            if T_start <= G * (tau - taui) + T_start <= T_end: 
                
                x = V * G * (tau - taui) + X_start

                t = G * (tau - taui) +  T_start
                break

            else:
                taui += (T_end-T_start)/G

        return [x, t]

    return traj



def tau(Ti, Xi):
    # turning points propertime
    Tau = [0]

    taui = 0
    for i, T in enumerate(Ti[:-1]):

        T_start = T 
        T_end = Ti[i + 1]

        X_start = Xi[i]
        X_end = Xi[i + 1]

        V = (X_end - X_start)/(T_end-T_start)

        G = 1/np.sqrt(1 - np.dot(V, V))

        taui += (T_end - T_start)/G

        Tau.append(taui)

    return Tau



def Sinc_peicewise(Xi, Ti):
    return 


def I(k, Omega, Xi, Ti):

    tau_turning_points = tau(Ti, Xi)

    tau_start = tau_turning_points[0]
    tau_end = tau_turning_points[-1]

    traj = peicewise(Xi, Ti)

    wave_vector = [-k, k]
    # already contracted with the metric

    def int_minus(tau):
        I_minus =  np.exp(1j * Omega * tau - 1j * np.dot(traj(tau), wave_vector))
        return I_minus
    
    def int_plus(tau):
        I_plus =  np.exp(1j * Omega * tau + 1j * np.dot(traj(tau), wave_vector))
        return I_plus

    return quad(int_minus, tau_start, tau_end)[0], quad(int_plus, tau_start, tau_end)[0]


k_grid =  np.linspace(-10, 10, 400)
Omega = 3

Xi = [0.1 * (-1)**i for i in range(20)]
Ti = [i for i in range(20)]

I_m = []
I_p = []

tau_turning_points = tau(Ti, Xi)
print(tau_turning_points)


tau_start = tau_turning_points[0]
tau_end = tau_turning_points[-1]

traj = peicewise(Xi, Ti)

t_test = np.linspace(tau_start, tau_end, 1000)

traj_test = np.array([traj(t) for t in t_test])

print(traj_test)


for k in k_grid:
    A, B = I(k, Omega, Xi, Ti)
    I_m.append(A)
    I_p.append(B)


fig, ax = plt.subplots(1, 2, figsize = [10, 5])
ax[0].plot(*traj_test.T, "-g")
ax[1].plot(k_grid, np.abs(I_m), "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p), "r", label="$I_+$")
ax[1].vlines(Omega, min(np.abs(I_m)), max(np.abs(I_m)), "r", label="$\Omega$", linestyle = "dashed")
ax[1].vlines(-Omega, min(np.abs(I_m)), max(np.abs(I_m)), "b", label="$-\Omega$", linestyle = "dashed")

ax[1].legend(frameon = False)
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")



plt.savefig("Wedge_traj.pdf")


# %%
