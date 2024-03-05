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

sp.print_latex(ARindler)

#%%
import sympy as sp 

g, X, rho, Q, r, eps  = sp.symbols(r"g X rho Q r e")

Phi = Q/r
Phi *= (1+g * X+g ** 2 * r ** 2/2)
Phi /= sp.sqrt(1 + g * X + g ** 2 * r ** 2/4)
phi = Phi/(1 + g * X)

Expansion_phi = sp.series(phi, g, 0, n=3)
Expansion_phi.simplify()

sp.print_latex(Expansion_phi)

Efield = sp.diff(phi.subs({r:sp.sqrt(X**2 + rho**2)}), X)
Expansion_phi2 = sp.series(Efield, X, sp.oo)
Expansion_phi2.simplify()
# %%
import sympy as sp 

g, X, rho, Q, r, eps  = sp.symbols(r"g X rho Q r e")

Phi = Q/r
Phi *= (1+g * X + g ** 2 * r ** 2/2)
Phi /= sp.sqrt(1 + g * X + g ** 2 * r ** 2/4)
Phi = Phi.subs({r:sp.sqrt(X**2 + rho**2)})


EX = sp.diff(Phi, X)
Erho = sp.diff(Phi, rho)
EX.simplify()
Erho.simplify()

Maxwell = sp.diff(EX * rho/(1+g*X), X) + sp.diff(Erho * rho/(1+g*X), rho)
Maxwell.simplify()


# %%
