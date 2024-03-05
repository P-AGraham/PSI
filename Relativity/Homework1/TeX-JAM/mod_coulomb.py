#%%
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
})

fig, ax = plt.subplots(1, 3, layout='constrained')

G = [0.4, 1 , 1.5]

mask = [0.05, 0.05, 0.05]

for i in range(3):
    g = G[i]
    delta = 1/(100*g)
    x = np.arange(-1/g-0.1, 1/g+0.1, delta)
    y = np.arange(-1/g-0.1, 1/g+0.1, delta)
    X, Y = np.meshgrid(x, y)

    ax[i].set_xticks([-1/g, 0, 1/g], ['$-g^{-1}$', '0', '$g^{-1}$'])

    



    Z1 = 1/np.sqrt(X**2 + Y**2)
    Z2 = (1 + g*X + g**2 * (X**2 + Y**2)/2)
    Z3 = np.sqrt(1 + g*X + g**2 * (X**2 + Y**2)/4)
    Z4 = 1/(1 + g * X)
    Z  = Z1 * Z2/(Z3 * Z4)

    interior = np.sqrt(X**2 + Y**2) < mask[i]
    Z[interior] = np.ma.masked

    nr, nc = Z.shape

    CS = ax[i].contourf(
    X, Y, Z,
    cmap=plt.cm.bone,
    levels=np.linspace(np.min(Z), np.max(Z)/g, 10)
    )

    ax[i].set_aspect('equal')

    if i == 0:
        ax[i].set_yticks([-1/g, 0, 1/g], ['$-g^{-1}$', '0', '$g^{-1}$'])
        x1, x2, y1, y2 = -0.3, 0.3, -0.3, 0.3  # subregion of the original image
        axins = ax[i].inset_axes(
        [0.1, 0.1, 0.3, 0.3],
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

        axins.contourf(
            X, Y, Z1 * Z2/(Z3 * Z4),
            cmap=plt.cm.bone,
            levels=np.linspace(np.min(Z), np.max(Z), 10)
            )
        
        ax[0].indicate_inset_zoom(axins, edgecolor="red")

    else:
        ax[i].set_yticks([], [])

    ax[i].axvline(x = -1/g, color = 'r', linestyle="dashed")

    ax[i].set_title('$g = {}$'.format(g))






ax[1].set_xlabel('$X$')
ax[0].set_ylabel(r'$\rho$')
ax[0].set_aspect('equal')


plt.savefig("Equi_Sketch.pdf", bbox_inches='tight')
# %%
