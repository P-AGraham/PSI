#%%
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
})


matplotlib.rcParams.update({'font.size': 22})

a, b = 1, 0.07

fig, ax = plt.subplots(1, 1, figsize=[5, 5], layout='constrained')

T = [3, 3.55, 3.9]

v = np.linspace(0.07001, 0.7, 1000)

for t in T: 
    P = t * np.exp(-a/(t * v))/(v-b)
    if t < 3.55:
        ax.plot(v, P, "b-", label="$T < T_c$")
    if t == 3.55:
        ax.plot(v, P, "r-", label="$T = T_c$")
    if t > 3.55:
        ax.plot(v, P, "k-", label="$T > T_c$")

ax.set_xlabel('$v$')
ax.set_ylabel('$P$')

ax.set_xticks([b], ['$b$'])
ax.set_yticks([0], ['$0$'])

ax.legend(frameon=False)
ax.set_ylim(3, 10)
#ax.set_aspect('equal')


plt.savefig("Isotherms.pdf", bbox_inches='tight')
# %%
