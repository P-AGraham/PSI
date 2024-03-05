#%%
import sympy as sp


Sx = sp.Matrix(([0,sp.sqrt(3),0,0],[sp.sqrt(3),0,2,0],[0,2,0,sp.sqrt(3)],[0,0,sp.sqrt(3),0]))
Sy = sp.Matrix(([0,-sp.sqrt(3),0,0],[sp.sqrt(3),0,-2,0],[0,2,0,-sp.sqrt(3)],[0,0,-sp.sqrt(3),0])) * sp.I
Sz = sp.Matrix(([3,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-3]))

T1 = 1/sp.sqrt(3)*( Sx + Sy + Sz)
T2 = 1/sp.sqrt(3)*( Sx - Sy - Sz)
T3 = 1/sp.sqrt(3)*(-Sx + Sy - Sz)
T4 = 1/sp.sqrt(3)*(-Sx - Sy + Sz)

T1.simplify()
T2.simplify()
T3.simplify()
T4.simplify()

T1 * T2 - T2 * T1

# https://en.wikipedia.org/wiki/Platonic_solid for cartesian coordinates
# %%
