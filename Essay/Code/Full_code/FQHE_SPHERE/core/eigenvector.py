from numpy import pi
import math, os, gc, copy
import numpy as np
from ctypes import *
from scipy.optimize import minimize, Bounds

from ..general import FQHE, abpath, mkdir
#from ..core import interaction
from ..tools.check import *
#from ..core.Hilbertspace.SLFermion import SLFermion_wrapper as cfl
#from ..core.Hamiltonian import csr, csr_fqhe



class eigenvector():
	def __init__(self, vec, eig, fqhe, qn, measured_qn={}):
		'''
			qn: quantum number
			example: qn = {"Lz":0}
					 qn = {"Lz":0, "Z_2":-1}
		'''
		if not isinstance(vec, np.ndarray):
			raise ValueError("unsupported type", type(vector))
		if not len(vec.shape)==1:
			raise ValueError("vector.shape = ", vector.shape)
		if not isinstance(fqhe, FQHE):
			raise ValueError(ins, "is not a FQHE isinstance!")
		self.__qn = qn.copy()
		self.__measured_qn = measured_qn.copy()
		self.__fqhe = fqhe.copy()
		self.__fqhe.clA()
		if not vec.flags['C_CONTIGUOUS']: 
			vec = np.ascontiguousarray(vec, dtype=vec.dtype)
		self.__vec = vec
		self.__eig = eig

	def __repr__(self):
		out = "The eigenvector in sector: "
		out += str((self.qn))
		out += "with E = " + str(self.eig) + ", dimension = " + str(self.shape[0]) + ".\n"
		out += "measured qn: " + str((self.measured_qn)) + "\n"
		return out

	def copy(self):
		return copy.deepcopy(self)

	def __lt__(self, obj):
		return self.eig < obj.eig

	@property
	def eig(self):
		return self.__eig

	@property
	def vec(self):
		return self.__vec

	@property
	def shape(self):
		return self.__vec.shape

	@property
	def qn(self):
		return self.__qn

	def get_qn(self, qn_name):
		if qn_name=="energy" and self.__qn.get(qn_name)==None:
			return self.__eig
		if qn_name=="eig":
			return self.__eig
		else:
			return self.__qn.get(qn_name)

	def get(self, qn_name):
		if self.get_qn( qn_name )!=None:
			return self.get_qn( qn_name )
		else:
			return self.get_measured_qn(qn_name)

	@property
	def measured_qn(self):
		return self.__measured_qn

	def get_measured_qn(self, qn_name):
		if self.__measured_qn.get(qn_name)!=None:
			return self.__measured_qn.get(qn_name)
		elif hasattr( self, qn_name):
			return getattr(self, qn_name)
		else:
			return None

	@property
	def fqhe(self):
		return self.__fqhe

	def add_measured_qn(self, qn_dict):
		self.__measured_qn.update(qn_dict.copy())
		return self.__measured_qn

	def dot(self, a):
		if isinstance(a, eigenvector):
			if self.qn != vec.qn:
				raise ValueError("self.qn = ", self.qn, "!= vec.qn = ", vec.qn)
			return np.vdot(self.vec, a.vec)
		return np.vdot(self.vec, a)

	@property
	def norm(self):
		return np.linalg.norm(self.vec)

	def __add__(self, vec):
		if not isinstance(vec, eigenvector):
			raise ValueError("vec.type = ", type(vec), "is not eigenvector!")
		if self.qn != vec.qn:
			raise ValueError("self.qn = ", self.qn, "!= vec.qn = ", vec.qn)
		return eigenvector(self.vec+vec.vec, "add", self.fqhe, self.qn)

	def __mul__(self, a):
		return eigenvector(self.vec*a, "mul", self.fqhe, self.qn)

	def __rmul__(self, a):
		return eigenvector(a*self.vec, "mul", self.fqhe, self.qn)

	def __truediv__(self, a):
		return eigenvector(self.vec/a, "div", self.fqhe, self.qn)

	@property
	def conj(self):
		return eigenvector(self.vec.conj(), "conj", self.fqhe, self.qn)

	@property
	def density(self):
		return self.fqhe.density(self)

	@property
	def sigma_x(self):
		return self.fqhe.sigma_x(self)

	@property
	def Z_2(self):
		return self.Z2

	@property
	def Z2(self):
		if self.qn.get("Z_2")!=None:
			return self.qn.get("Z_2")
		elif self.measured_qn.get("Z_2")!=None:
			return self.measured_qn.get("Z_2")
		resu = self.fqhe.Z2(self)
		self.add_measured_qn({"Z_2":resu})
		return resu

	@property
	def PH(self):
		if self.qn.get("PH")!=None:
			return self.qn.get("PH")
		elif self.measured_qn.get("PH")!=None:
			return self.measured_qn.get("PH")
		resu = self.fqhe.PH(self)
		self.add_measured_qn({"PH":resu})
		return resu