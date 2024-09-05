#%%
import numpy as np 
from scipy.integrate import quad 
import matplotlib.pyplot as plt

def make_traj(v, T):

    g = 1/np.sqrt(1-v**2)

    def traj(tau):

        if tau <= -T/2:
            return 0, tau
        if tau > T/2:
            return 0, tau
        if -T/2 < tau <= 0:
            return g*v*(T/2 + tau), g*tau
        if 0 < tau <= T/2:
            return g*v*(T/2 - tau), g*tau
    return traj


traj = make_traj(0.3, 1)

Tau = np.linspace(-2, 2)
XT = np.array([traj(tau) for tau in Tau])
X = XT[:, 0]
T = XT[:, 1]

plt.plot(X, T, ".-")

#%%

def I_direct(k, Omega, traj, start, end):
    kv = np.abs(k)

    def int_minus(tau):
        I_minus =  np.exp(1j * Omega * tau - 1j * traj(tau)[1] * kv) * np.sin(traj(tau)[0] * k)
        return I_minus
    
    def int_plus(tau):
        I_plus =  np.exp(1j * Omega * tau + 1j * traj(tau)[1] * kv) * np.sin(-traj(tau)[0] * k)
        return I_plus
    

    return quad(int_minus, start, end)[0], quad(int_plus, start, end)[0]

T = 110

k_grid =  np.linspace(20, 200, 300)
Omega = 1
I_m, I_p  = [], []

traj = make_traj(0.1, 1)

t_test = np.linspace(-2*T, 2*T, 100)
traj_test = np.array([traj(t) for t in t_test])

for k in k_grid:
    A, B = I_direct(k, Omega, traj, -2*T, 2*T)
    I_m.append(A)
    I_p.append(B)



fig, ax = plt.subplots(1, 2, figsize = [10, 5])
ax[0].plot(*traj_test.T, "-g")
ax[1].plot(k_grid, np.abs(I_m)**2, "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p)**2, "r--", label="$I_+$")
#ax[1].set_xlim(75, 76)
#ax[1].set_ylim(0, 1)
ax[0].legend(frameon = False, fontsize="15")
ax[1].legend(frameon = False, fontsize="15")
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")
# %%
