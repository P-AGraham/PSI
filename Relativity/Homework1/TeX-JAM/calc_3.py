#%%
import sympy as sp 

x, t, rho, phi, X, T, g, Q = sp.symbols(r"x t \rho \varphi X T g Q")

L = 1/g
xi = sp.sqrt((L**2 + t**2 - rho**2 - x**2)**2/4 + L**2 * rho**2)
delta = rho**2 + x**2 + L**2 - t**2

xQ = (x * delta - 2 * t * xi)/(2 * x ** 2 - 2 * t ** 2)
tQ = (t * delta - 2 * x * xi)/(2 * x ** 2 - 2 * t ** 2)

xc = (X+1/g) * sp.cosh(g * T)
tc = (X+1/g) * sp.sinh(g * T)

A = sp.Matrix([-xQ * Q/xi, tQ * Q/xi, 0, 0])
A.simplify()
A = A.subs({x : xc, t : tc})
A.simplify()

J = sp.Matrix([tc, xc, rho, phi]).jacobian([T, X, rho, phi]).T

ARindler = J * A
ARindler.simplify()
ARindler

# %%
# %%
