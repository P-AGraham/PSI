#%%
import matplotlib.pyplot as plt
import numpy as np
plt.style.use(["tableau-colorblind10",
                "C:/Users/pgraham1/Documents/GitHub/PSI/Quantum_Matter/Homework1/Essay.mplstyle"])


plt.rcParams.update({
    "text.usetex": True,
})

fig, ax = plt.subplots(1, 1, layout='constrained', figsize=[7, 4])

@np.vectorize
def ev_Z(gJ):
    if gJ/2 < 1:
        return np.sqrt(1-(gJ/2)**2)
    else: 
        return 0

gJ_grid =  np.linspace(0, 4, 500)

ax.plot(gJ_grid, ev_Z(gJ_grid), "b-", linewidth=2)

ax.set_xlabel('$g/J$')
ax.set_ylabel(r'$\langle Z_i \rangle$')
#ax.set_aspect('equal')

plt.axvline(x=2, ymin=0, color='r', linestyle='--', label="Phase transition")
ax.text(0.15, 0.35, 'Ferromagnet', fontsize=16)
ax.text(0.15+2.5, 0.35, 'Paramagnet', fontsize=16)

ax.legend(fontsize=13, frameon=False)

plt.savefig("Z_average.pdf", bbox_inches='tight')

# %%
