#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import math
from functools import reduce
import scipy.linalg as sla
import scipy.sparse as sp
import time


#%%

Sz = sp.dia_matrix([[1.0,0.0],[0.0,-1.0]])
Sx = sp.csc_matrix([[0.0,1.0],[1.0,0.0]])
Sy = sp.csc_matrix([[0.0,-1j],[1j,0.0]])

def Id(n): #defining the identity matrix
	return sp.eye(2**n)

def Szi(i,n): #defining the Sz operator at site i given n spins
	A = Id(i)
	B = Id(n-i-1)
	D = reduce(sp.kron, [A,Sz,B])
	return D

def Sxi(i,n): #defining the Sx operator at site i given n spins
	A = Id(i)
	B = Id(n-i-1)
	D = reduce(sp.kron, [A,Sx,B])
	return D

def interaction(c,i,j,n): #defining the SzSz operator at sites i, j given n spins
	if (i == j):
		print("i and j must be distinct!")
	elif (j<i):
		return interaction(c,j,i,n)
	else:
		A = Id(i)
		B = Id(j-i-1)
		C = Id(n-j-1)
		D = reduce(sp.kron, [A,Sz,B,Sz,C])
	return c*D


def NNterm(n): 
    '''Using Periodic BC'''
    for i in range(n):
        if(i == 0):
            H = interaction(1,i,i+1,n)
        else:
            H = H + interaction(1,i,(i+1)%n,n)
    return H

def TrField(n): #defining the transverse magnetic field term
	for i in range(n):
		if (i == 0):
			F = Sxi(i,n)
		else:
			F = F + Sxi(i,n)
	return F


def sparse_diag(n, k, g = 1):
	"""
	a function that provides the sparse diagonalization of the 1D TFIM model given a system size n and the transverse field g
	returns: (energies, eigenvectors)
	
	g is set to 1 because we study the CFT at the quantum critical point
	
	k is the number of low energy eigenstates extracted
	"""

	trfield = TrField(n)
	interterm = NNterm(n)

	H = -sp.csc_matrix(interterm + g*trfield) #sparse representation of the Hamiltonian
	
	val, vecs = sp.linalg.eigsh(H,k=k,which='SA')
	return val, vecs, H

def Hmode(n, m, g=1):
    """
    A function calculating the Fourier mode m of H used in identification of primaries for the TFIM of size n
    """
    N = n
    Hm = - np.exp(1j * (0 + 1/2) * m * 2 * np.pi/N) * interaction(1,0,1,n)
    Hm = Hm - g * np.exp(1j * 0 * m * 2 * np.pi/N) * Sxi(0,n)

    for i in range(1, n):
        Hm = Hm - np.exp(1j * (i + 1/2) * m * 2 * np.pi/ N) * interaction(1,i,(i+1)%n,n) 
        Hm = Hm - g * np.exp(1j * i * m * 2 * np.pi/N) * Sxi(i,n)

    Hm = N * Hm/(2 * np.pi)

    return Hm


def rescale_shift(n, k, g=1): 
	"""
	a function applying an afine transformation to the TFIM Hamiltonian to bring its spectrum to match approximatly the CFT Hamiltonian spectrum
	"""

	val, vecs, H = sparse_diag(n, k, g=g)

	# Shifting the ground state to 0 (corresponds to the ientity operator which has scaling dimension 0)
 
	val_shifted = val - val[0] 

	# Finding the rescaling factor making the spectrum consistent with CFT scaling dimensions at low energy 
 
	O2 = (Hmode(n, 2, g=g) + Hmode(n, -2, g=g))/2

	ground_state = vecs[:, 0]

	Delta_I2 = 2

	new_H = H - val[0] * sp.csc_matrix(Id(n))

	num = O2 @ ground_state
	den = new_H @ num 

	a = np.linalg.norm(num) * Delta_I2/np.linalg.norm(den) 
	
	# Rescaling the spectrum 
 
	rescaled_shifted_val = val_shifted * a

	b = - a * val[0]

	new_H = a * new_H
	return rescaled_shifted_val, vecs, a, b, new_H

def primary_test(n, k, eps_max = 1e-12, g = 1):
	"""
	a function testing all k lowest energy eigenstates of the TFIM of size n to find primary candidates consistent with the threshold eps_max
	"""
	val, vecs, a, b, _ = rescale_shift(n, k, g=g)

	print(a, b)

	O1 = (Hmode(n, 1, g=g) + Hmode(n, -1, g=g))/2
	O2 = (Hmode(n, 2, g=g) + Hmode(n, -2, g=g))/2

	test = []
	for i in range(k):

		state1_0 = O1 @ vecs[:, i]
		state2_0 = O2 @ vecs[:, i]

		state1 = np.zeros(state1_0.shape, dtype="complex128")
		state2 = np.zeros(state2_0.shape, dtype="complex128")

		if i == 0: 
			test.append(True)

		else:
			for state_idx in range(i):
				v = vecs[:, state_idx]

				state1 += np.dot(np.conj(v),  state1_0) * v
				state2 += np.dot(np.conj(v), state2_0) * v
				

			eps = np.linalg.norm(state1)
			eps = eps + np.linalg.norm(state2)

			test.append(eps)

	return np.round(val, 4), np.round(np.log(test), 2)


#%%

# Shifted and rescaled spectrum
spectrum, eps = primary_test(16, 22)

print(spectrum, eps)

spectrum2, eps = primary_test(18, 3)

# Critical exponents calculations

D = 2 # spacetime dimension

Delta_sigma = spectrum[1]
Delta_eps = spectrum[2]

Delta_sigma2 = spectrum2[1]
Delta_eps2 = spectrum2[2]

## alpha 

alpha = 2-D/(D-Delta_eps)
alpha2 = 2-D/(D-Delta_eps)

## beta

beta = Delta_sigma/(D - Delta_eps)
beta2 = Delta_sigma2/(D - Delta_eps2)

## gamma

gamma = (D - Delta_sigma)/(D - Delta_eps)
gamma2 = (D - Delta_sigma2)/(D - Delta_eps2) 

## delta

delta = (D - Delta_sigma)/Delta_sigma
delta2 = (D - Delta_sigma2)/Delta_sigma2

## eta

eta = 2 * Delta_sigma - D + 2
eta2 = 2 * Delta_sigma2 - D + 2

## nu

nu = 1/(D - Delta_eps)
nu2 = 1/(D - Delta_eps2)
 

exp_exact = [0, 1/8, 7/4, 15, 1/4, 1]

exp_1 = [alpha, beta, gamma, delta, eta, nu]
exp_2 = [alpha2, beta2, gamma2, delta2, eta2, nu2]

print(exp_exact, exp_1, exp_2)


# %%
# Numerical identification of descendant operators
val, vecs, _, _, H = rescale_shift(16, 22, g = 1)

ket_I = vecs[:, 0]
ket_sigma = vecs[:, 1]
ket_eps = vecs[:, 2]

O1 = (Hmode(16, 1) + Hmode(16, -1))/2
O2 = (Hmode(16, 2) + Hmode(16, -2))/2


# sigma descendants

den_sigma1 = O1 @ ket_sigma
den_sigma2 = O2 @ ket_sigma

num_sigma1 = H @ den_sigma1
num_sigma2 = H @ den_sigma2

Delta_sigma1 = np.linalg.norm(num_sigma1)/np.linalg.norm(den_sigma1)
Delta_sigma2 = np.linalg.norm(num_sigma2)/np.linalg.norm(den_sigma2)

Delta_sigma1_exact =  val[1] + 1
Delta_sigma2_exact =  val[1] + 2

# epsilon descendants

den_eps1 = O1 @ ket_eps
den_eps2 = O2 @ ket_eps

num_eps1 = H @ den_eps1
num_eps2 = H @ den_eps2

Delta_eps1 = np.linalg.norm(num_eps1)/np.linalg.norm(den_eps1)
Delta_eps2 = np.linalg.norm(num_eps2)/np.linalg.norm(den_eps2)
#%%
Delta_eps1_exact =  val[2] + 1
Delta_eps2_exact =  val[2] + 2


# reuslts 

c1 = [Delta_sigma1, Delta_sigma2]
c2 = [Delta_sigma1_exact, Delta_sigma2_exact]
c3 = [Delta_eps1, Delta_eps2]
c4 = [Delta_eps1_exact, Delta_eps2_exact]


# %%
from scipy.optimize import curve_fit

def c_model(x, a, b): 
    return a*x + b

# Central charge culculation through finite size scaling 
plt.style.use(["tableau-colorblind10",
                "C:/Users/pgraham1/Documents/GitHub/PSI/Quantum_Matter/Homework1/Essay.mplstyle"])

plt.rcParams.update({
    "text.usetex": True,
})

fig, ax = plt.subplots(1, 1, layout='constrained', figsize=[7, 4])

L_list = np.arange(7, 17, 1)
L_space =  np.linspace(7, 17)

c_list = []
for L in L_list:
	val, vecs, a, b, H = rescale_shift(L, 1)
	ground_state = vecs[:, 0]

	H2 = Hmode(L, -2)
	
	c = np.linalg.norm(H2 @ ground_state)**2/2
	c_list.append(c)


popt, pcov = curve_fit(c_model, 1/L_list**2, c_list) 

ax.plot(1/L_space**2, c_model(1/L_space**2, *popt),
		 label=r"$c_L = a L^{-2} + c_{\infty}$"
		 )

ax.plot([], [], "w.",
		 label=r"$c_\infty \approx $" + "${:.5f}$".format(popt[1])
		 )

ax.plot(1/L_list**2, c_list, ".")


ax.set_xlabel("$L^{-2}$")
ax.set_ylabel("$c_L$")
ax.legend(frameon=False, fontsize=16)

plt.savefig("central_charge_fit.pdf")



# %%
