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

    #wave_vector = [k, -k]

    # already contracted with the metric

    kv = np.abs(k)

    def int_minus(tau):
        I_minus =  np.exp(1j * Omega * tau - 1j * traj(tau)[1] * kv) * np.sin(traj(tau)[0] * k)
        return I_minus
    
    def int_plus(tau):
        I_plus =  np.exp(1j * Omega * tau + 1j * traj(tau)[1] * kv) * np.sin(-traj(tau)[0] * k)
        return I_plus
    

    #def int_minus(tau):
    #    I_minus = np.exp(1j * Omega * tau - 1j * traj(tau)[0] * kv + 1j * traj(tau)[1] * k)
    #    return I_minus
    #
    #def int_plus(tau):
    #    I_plus = np.exp(1j * Omega * tau + 1j * traj(tau)[0] * kv - 1j * traj(tau)[1] * k)
    #    return I_plus
    

    I_minus_int = 0
    I_plus_int = 0
    for i, tp in enumerate(tau_turning_points[:-1]):
        start = tp 
        end = tau_turning_points[i+1]
        I_minus_int += quad(int_minus, start, end)[0]
        I_plus_int += quad(int_plus, start, end)[0]
        
    return I_minus_int, I_plus_int 




def alpha_dot_wrap(T, v0, v, k):

    a = (v - v0)/T

    def alpha_dot(tau):
        if tau < 0:
            return  k * v0
        elif 0 <= tau < T:
            return  k * (v0 + a * tau)
        elif tau >= T:
            return  k * v 

    return alpha_dot


def I_alpha(Omega, start, end, k = 1, T = 5.105, v = 1.801, v0 = 1.03):

    alpha_dot = alpha_dot_wrap(T, v0, v, k)

    def alpha(tau):
        return quad(alpha_dot, 0, tau)[0]

    def int_minus(tau):
        I_minus = np.exp(1j * Omega * tau - 1j * alpha(tau))
        return I_minus
    
    def int_plus(tau):
        I_plus = np.exp(1j * Omega * tau + 1j * alpha(tau))
        return I_plus
        

    return quad(int_minus, start, end)[0], quad(int_plus, start, end)[0]


#Omega_grid =  np.linspace(0, 2, 200)
#I_m, I_p  = [], []
#start = 0
#end = 5.105
#
#for Omega in Omega_grid:
#    A, B = I_alpha(Omega, start, end)
#    I_m.append(A)
#    I_p.append(B)
#
#
#fig, ax = plt.subplots(1, 1, figsize = [10, 5])
#ax.plot(Omega_grid, np.log(np.abs(I_m)), "b", label="$I_-$")
#ax.plot(Omega_grid, np.log(np.abs(I_p)), "r", label="$I_+$")
#
#ax.legend(frameon = False)
#ax.set_title("Trajectory in $k$-plane")
#ax.set_title("$I_\pm(k)$")
#ax.set_xlabel("$\Omega$")

def circle_traj(omega, r):
    def traj(tau):
        gamma = 1/np.sqrt(1-r**2 * omega**2)
        x = gamma * r * np.cos(omega * tau)
        t = gamma * tau 
        return x, t
    return traj



def I_direct(k, Omega, traj, start, end):
    kv = np.abs(k)

    def int_minus(tau):
        I_minus =  np.exp(1j * Omega * tau - 1j * traj(tau)[1] * kv) * np.sin(traj(tau)[0] * k)
        return I_minus
    
    def int_plus(tau):
        I_plus =  np.exp(1j * Omega * tau + 1j * traj(tau)[1] * kv) * np.sin(-traj(tau)[0] * k)
        return I_plus
    

    return quad(int_minus, start, end)[0], quad(int_plus, start, end)[0]

# %%
