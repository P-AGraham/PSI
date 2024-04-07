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

	num = O2 @ ground_state
	den = (H - val[0] * sp.csc_matrix(Id(n))) @ num 

	a = np.linalg.norm(num) * Delta_I2/np.linalg.norm(den) 
	
	# Rescaling the spectrum 
 
	rescaled_shifted_val = val * a

	b = - a * val[0]

	return rescaled_shifted_val, vecs, a, b

def primary_test(n, k, eps_max = 1e-12, g = 1):
	"""
	a function testing all k lowest energy eigenstates of the TFIM of size n to find primary candidates consistent with the threshold eps_max
	"""
	val, vecs, a, b = rescale_shift(n, k, g=g)

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

			test.append(eps < eps_max)

	return dict(zip((val-val[0]), test))


#%%
primary_test(16, 22)


# %%

# %%
