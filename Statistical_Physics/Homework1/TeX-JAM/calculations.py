#%%
import sympy as sp 

T, a, b, v, k = sp.symbols("T a b v k")

ex = sp.exp(-a/(k * T * v))/(v-b)
P = k * T * ex

del1P = sp.diff(P, v, 1).simplify()
Expr = del1P.simplify() 

sp.print_latex(Expr)

# %%
import sympy as sp 

vr, Tr, x, y, t, c = sp.symbols("vr Tr x y t c")

P = sp.exp(2-2/(Tr * vr))/(2 * vr - 1) * Tr

P = P.subs({Tr:1-t}).series(t, 0, 2).removeO()

Pl = P.subs({vr:1-x}).series(x, 0, 4) 
Pg = P.subs({vr:1+y}).series(y, 0, 4) 
Pg = Pg.subs({y : x + c * x ** 2})

Eq = (Pg-Pl).simplify().removeO()

sp.print_latex(Eq)

X = sp.solve(Eq, x)[1].series(t, 0, 2)

sp.print_latex(X)

# %%
