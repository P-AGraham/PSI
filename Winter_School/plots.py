#%% 
from Chirp_v import *

#%%
# Jigsaw trajectory

Xi = [0 for i in range(4)]
Ti = [i for i in range(4)]


print(Ti)
Omega = 5
k_grid =  np.linspace(0, 10, 400)
I_m, I_p  = [], []

tau_turning_points = tau(Ti, Xi)
tau_start = tau_turning_points[0]
tau_end = tau_turning_points[-1]
traj = peicewise(Xi, Ti)
t_test = np.linspace(tau_start, tau_end, 1000)
traj_test = np.array([traj(t) for t in t_test])

for k in k_grid:
    A, B = I(k, Omega, Xi, Ti)
    I_m.append(B)
    I_p.append(A)


fig, ax = plt.subplots(1, 2, figsize = [10, 5])
ax[0].plot(*traj_test.T, "-g")
ax[1].plot(k_grid, np.abs(I_m)**2, "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p)**2, "r", label="$I_+$")

ax[1].legend(frameon = False)
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")

#ax[1].set_ylim(0, 0.2)

#ax[1].set_ylim(0, 0.5)
#
#ax[1].set_xlim(Omega-0.5, Omega+0.5)

plt.savefig("Figures/jigsaw_traj.pdf")
#%%





# %%
# Inertial trajectory

Xi = [(-1)**i*0.01 for i in range(11)]
Ti = [i for i in range(11)]

k_grid =  np.linspace(-15, 15, 400)
Omega = 8
I_m, I_p  = [], []

tau_turning_points = tau(Ti, Xi)
tau_start = tau_turning_points[0]
tau_end = tau_turning_points[-1]
traj = peicewise(Xi, Ti)
t_test = np.linspace(tau_start, tau_end, 1000)
traj_test = np.array([traj(t) for t in t_test])

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

plt.savefig("Figures/Inertial_traj.pdf")

# %%

# %%
# Hyperbolic
delta_t = 0.1
Xi = [np.cosh(i * delta_t) for i in range(10)]
Ti = [np.sinh(i * delta_t)  for i in range(10)]

Xi = Xi + Xi[::-1]
Ti = Ti + [np.sinh(9.5 * delta_t) + np.sinh(i * delta_t)  for i in range(10)]


k_grid =  np.linspace(-30, 30, 400)
Omega = 10
I_m, I_p  = [], []

tau_turning_points = tau(Ti, Xi)
tau_start = tau_turning_points[0]
tau_end = tau_turning_points[-1]
traj = peicewise(Xi, Ti)
t_test = np.linspace(tau_start, tau_end, 1000)
traj_test = np.array([traj(t) for t in t_test])

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

plt.savefig("Figures/Hyperbolic_traj.pdf")
# %%
# Broken 
delta_t = 0.1
Xi = [-1, -2, -2, -2]
Ti = [-2, -1, 1, 2]

k_grid =  np.linspace(-20, 20, 400)
Omega = 10
I_m, I_p  = [], []

tau_turning_points = tau(Ti, Xi)
tau_start = tau_turning_points[0]
tau_end = tau_turning_points[-1]
traj = peicewise(Xi, Ti)
t_test = np.linspace(tau_start, tau_end, 1000)
traj_test = np.array([traj(t) for t in t_test])

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

#ax[1].set_ylim(0, 5)

plt.savefig("Figures/Circular_traj.pdf")
# %%
