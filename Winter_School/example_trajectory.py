#%%
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import quad 
plt.style.use(["tableau-colorblind10", "C:/Users/pgraham1/Documents/GitHub/PSI/Winter_School/presentation.mplstyle"])



def alpha_dot_wrap(T1, T2, s0, s1, s2):

    def alpha_dot(tau):
        if tau < 0:
            return  s0
        elif 0 <= tau < T1:
            return s0 + (s1-s0)/T1 * tau
        elif 0 <= tau < T2:
            return s1 + (s2-s1)/T1 * tau
        elif tau >= T2:
            return  s2
        
    return alpha_dot

alpha_dot = alpha_dot_wrap(0.7, 1.5, 1, 1.1, 1.5)

def xmu_dot_wrap(alpha_dot, k):

    def t_dot(tau):
        return (k/alpha_dot(tau) + alpha_dot(tau)/k)/2
    
    def x_dot(tau):
        return (k/alpha_dot(tau) - alpha_dot(tau)/k)/2
    
    return t_dot, x_dot


tau0 = -0.4
tau_grid = np.linspace(tau0, 1.5)

t_dot, x_dot = xmu_dot_wrap(alpha_dot, 1)

t = []
x = []

for taustep in tau_grid:
    t.append(quad(t_dot, tau0, taustep)[0])
    x.append(quad(x_dot, tau0, taustep)[0])

fig, ax = plt.subplots(1, 1, figsize = [4, 8])
ax.plot(x, t, "-g")

ax.set_title("Trajectory in $k$-plane")
ax.set_xlabel("$x$")
ax.set_ylabel("$ct$")

ax.set_xlim(-0.05, 0.002)
ax.set_ylim(0, 1.2)



# %%
