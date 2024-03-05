#%%
from Chirp_v import *
plt.style.use(["tableau-colorblind10", "C:/Users/pgraham1/Documents/GitHub/PSI/Winter_School/presentation.mplstyle"])
#%%
# Jigsaw trajectory

v = 0.2973
Xi = [v*(-1)**(i) for i in range(9)]
Ti = [i for i in range(9)]


print(Ti)
Omega = 3
k_grid =  np.linspace(21.8, 22.5, 200)
k_grid_zoom =  np.linspace(7, 10, 200)
I_m, I_p  = [], []
I_m_zoom, I_p_zoom = [], []

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


for k in k_grid_zoom:
    A, B = I(k, Omega, Xi, Ti)
    I_m_zoom.append(A)
    I_p_zoom.append(B)

L = 0.34

max_n = int(np.round((max(k_grid) * L/np.pi)))
min_n = int(np.round((min(k_grid) * L/np.pi)))+1
cavity = [np.pi*n/L for n in range(min_n, max_n)]


fig, ax = plt.subplots(1, 2, figsize = [12, 5])
ax[0].plot(*traj_test.T, "-g", label="$v=${}$c$".format(v*2))
ax[1].plot(k_grid, np.abs(I_m)**2, "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p)**2, "r--", label="$I_+$")

ax[1].plot([22.07], [min(np.abs(I_m)**2)], "kX", label="Cavity mode")

ax[1].legend(frameon = False, fontsize="15")
ax[0].legend(frameon = False, fontsize="15")
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$|I_\pm(k)|^2$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")

ax[1].set_yscale('log')

#ax[1].set_ylim(0, 0.2)

#ax[1].set_ylim(0, 0.5)
#
#ax[1].set_xlim(Omega-0.5, Omega+0.5)

plt.savefig("Figures/jigsaw_traj_test_look_like_circular_shift_3.pdf")
# %%

Xi = [0.299*(-1)**(i) for i in range(9)]
Ti = [i for i in range(9)]


print(Ti)
Omega = 3
k_grid =  np.linspace(22, 23, 500)
k_grid_zoom =  np.linspace(7, 10, 500)
I_m, I_p  = [], []
I_m_zoom, I_p_zoom = [], []

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


for k in k_grid_zoom:
    A, B = I(k, Omega, Xi, Ti)
    I_m_zoom.append(A)
    I_p_zoom.append(B)

L = 0.34

max_n = int(np.round((max(k_grid) * L/np.pi)))
min_n = int(np.round((min(k_grid) * L/np.pi)))+1
cavity = [np.pi*n/L for n in range(min_n, max_n)]


fig, ax = plt.subplots(1, 2, figsize = [12, 5])
ax[0].plot(*traj_test.T, "-g", label="$v=${}$c$".format(0.299*2))
ax[1].plot(k_grid, np.abs(I_m), "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p), "r--", label="$I_+$")

ax[1].set_yscale('log')

#ax[1].plot(cavity, [0]*(max_n-min_n), "mX", label="Cavity modes")

ax[1].legend(frameon = False, fontsize="15")
ax[0].legend(frameon = False, fontsize="15")
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")


#ax[1].set_ylim(0, 0.0002)

#ax[1].set_ylim(0, 0.5)
#
#ax[1].set_xlim(Omega-0.5, Omega+0.5)

plt.savefig("Figures/jigsaw_traj_quartic_good9cycles_2.pdf")

# %%
Xi = [0.31*(-1)**(i) for i in range(3)]
Ti = [i for i in range(3)]


print(Ti)
Omega = 3
k_grid =  np.linspace(22, 23, 500)
k_grid_zoom =  np.linspace(7, 10, 500)
I_m, I_p  = [], []
I_m_zoom, I_p_zoom = [], []

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


for k in k_grid_zoom:
    A, B = I(k, Omega, Xi, Ti)
    I_m_zoom.append(A)
    I_p_zoom.append(B)

L = 0.34

max_n = int(np.round((max(k_grid) * L/np.pi)))
min_n = int(np.round((min(k_grid) * L/np.pi)))+1
cavity = [np.pi*n/L for n in range(min_n, max_n)]


fig, ax = plt.subplots(1, 2, figsize = [12, 5])
ax[0].plot(*traj_test.T, "-g")
ax[1].plot(k_grid, np.abs(I_m), "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p), "r--", label="$I_+$")

ax[1].set_yscale('log')

#ax[1].plot(cavity, [0]*(max_n-min_n), "mX", label="Cavity modes")

ax[1].legend(frameon = False, fontsize="15")
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")

#ax[1].set_ylim(0, 0.0002)

#ax[1].set_ylim(0, 0.5)
#
#ax[1].set_xlim(Omega-0.5, Omega+0.5)

plt.savefig("Figures/jigsaw_traj_quartic_good1cycles.pdf")
# %%
# Circular

#r = 0.15
r = 0.14579 # Robust
T = 2

omega= 2*np.pi/T

traj = circle_traj(omega, r)


k_grid =  np.linspace(19.5, 19.8, 500)
Omega = 1
I_m, I_p  = [], []

t_test = np.linspace(0*T, 5*T, 1000)
traj_test = np.array([traj(t) for t in t_test])

for k in k_grid:
    A, B = I_direct(k, Omega, traj, 0*T, 5*T)
    I_m.append(A)
    I_p.append(B)



fig, ax = plt.subplots(1, 2, figsize = [10, 5])
ax[0].plot(*traj_test.T, "-g", label="$\omega r = {}c$".format(np.round(omega*r, 3)))
ax[1].plot(k_grid, np.abs(I_m)**2, "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p)**2, "r", label="$I_+$")

ax[0].legend(frameon = False, fontsize="15")
ax[1].legend(frameon = False, fontsize="15")
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")

ax[1].set_yscale('log')

#ax[1].set_ylim(0, 0.0002)

#ax[1].set_ylim(0, 0.5)
#
#ax[1].set_xlim(Omega-0.5, Omega+0.5)

plt.savefig("Figures/circular_traj_quartic_tests.pdf")
# %%
# mode alignement 

r = 0.05
T = 2

omega= 2*np.pi/T

traj = circle_traj(omega, r)


k_grid =  np.linspace(0, 20, 300)
Omega = 1
I_m, I_p  = [], []

t_test = np.linspace(0*T, 1*T, 1000)
traj_test = np.array([traj(t) for t in t_test])

for k in k_grid:
    A, B = I_direct(k, Omega, traj, 0*T, 1*T)
    I_m.append(A)
    I_p.append(B)



fig, ax = plt.subplots(1, 2, figsize = [10, 5])
ax[0].plot(*traj_test.T, "-g", label="$\omega r = {}c$".format(np.round(omega*r, 3)))
ax[1].plot(k_grid, np.abs(I_m), "b", label="$I_-$")
ax[1].plot(k_grid, np.abs(I_p), "r--", label="$I_+$")

ax[0].legend(frameon = False, fontsize="15")
ax[1].legend(frameon = False, fontsize="15")
ax[0].set_title("Trajectory in $k$-plane")
ax[1].set_title("$I_\pm(k)$")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$ct$")
ax[1].set_xlabel("$k$")

#x1, x2, y1, y2 = 5, 8, 0, 0.01  # subregion of the original image
#axins = ax[1].inset_axes(
#    [0.5, 0.2, 0.47, 0.47],
#    xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
#axins.plot(k_grid, np.abs(I_m)**2, "b")
#axins.plot(k_grid, np.abs(I_p)**2, "r--")
#axins.plot(cavity, [0]*(max_n-min_n), "mX")

#ax[1].indicate_inset_zoom(axins, edgecolor="black")

#ax[1].set_yscale('log')

#ax[1].set_ylim(0, 0.0002)

#ax[1].set_ylim(0, 0.5)
#
#ax[1].set_xlim(Omega-0.5, Omega+0.5)

plt.savefig("Figures/alignement_example.pdf")
# %%
