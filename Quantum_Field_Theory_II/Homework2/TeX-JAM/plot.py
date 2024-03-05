#%%
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
})

matplotlib.rcParams.update({'font.size': 22})

fig, ax = plt.subplots(1, 1, figsize = [7, 5],  layout='constrained')


gr = np.linspace(0, 0.5)



ax.set_xlabel('$g_R(\mu)$')

ax.plot(gr, 3 * gr**2/(4 * np.pi)**2, "r-")
ax.set_ylabel(r'$\beta_g(g_R(\mu))/\hbar$')
#ax.set_aspect('equal')


plt.savefig("beta.pdf", bbox_inches='tight')
# %%
