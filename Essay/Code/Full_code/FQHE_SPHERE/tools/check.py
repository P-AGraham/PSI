
#from ..fermion.singlelayer import FermionSL

def ispositiveint(n):
	if isinstance(n, int) and n >= 0:
		return True
	else:
		return False

def isFermionSL(ins):
	if isinstance(ins, FermionSL):
		return True
	else:
		return False

def isFermionSL(ins):
	if isinstance(ins, FermionTL):
		return True
	else:
		return False

def gcd_lcm(x, y):
	if ispositiveint(x) and ispositiveint(y):
		a, b = x, y
		while a != b:
			if a > b:
				a -= b
			else:
				b -= a
		return a, int(x*y/a) #(gcd, lcm)
	else:
		raise ValueError("(x, y) must be positive integers!")