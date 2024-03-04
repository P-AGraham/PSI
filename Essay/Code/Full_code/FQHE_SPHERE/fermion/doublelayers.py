import math, os, gc, copy
import numpy as np
from ctypes import *
from numpy import pi
from scipy.special import sph_harm

import time
from ..tools.check import ispositiveint, gcd_lcm
from ..general import FQHE, abpath, mkdir
from ..core import interaction
from ..core.Hamiltonian import csr, csr_fqhe, csr_H_L2
from ..core.eigenvector import eigenvector as eigvec
#from ..core.eigenvector import eigenvector_k1
from ..core.Hilbertspace.QHFMFermion import QHFMFermion_wrapper as cfl
#from ..core.Hilbertspace.SLFermion import SLFermion_wrapper as cfl_inter

indtype = np.uint64
class FermionDL(FQHE):
	def __init__(self, n_e, n_o, num_thread=8):
		super().__init__(n_e, n_o, num_thread)

'''
class FermionDL(FQHE):
	def __init__(self, n_phi=16, n_up=4, n_down=4, tau=(0.0, 1.0), num_thread=6):
		super().__init__(n_up+n_down, tau, num_thread)
		self.__level = 0
		if ispositiveint(n_phi):
			self.__n_phi = n_phi
		else:
			raise ValueError("The orbital number n must be a positive integer!")
		if ispositiveint(n_up) and ispositiveint(n_down):
			self.__n_up = n_up
			self.__n_down = n_down
		else:
			raise ValueError("The particle number n_up or n_down must be a positive integer!")

		self.__area = 2 * pi * self.__n_phi
		self.__lx = math.sqrt(self.__area/self.tau[1])
		self.__ly = self.__lx*math.sqrt(self.tau[0]**2+self.tau[1]**2)
		self.__theta = np.angle(self.tau[0]+1j*self.tau[1])
		cfl.init(self.__n_phi, self.__n_up, self.__n_down, self.num_thread)
		self.__A1_2body = 0
		self.__A2_2body = 0
		self.__A12_2body= 0
		self.__pqgcd = gcd_lcm(n_phi, self.n)[0]
		self.__q = n_phi//self.__pqgcd
		self.__k1, self.__k2 = -1, -1

	def __repr__(self):
		return str(self.geometry)

	def deformation_tau(self, tau):
		super().deformation_tau(tau)
		self.__lx = math.sqrt(self.__area/self.tau[1])
		self.__ly = self.__lx*math.sqrt(self.tau[0]**2+self.tau[1]**2)
		self.__theta = np.angle(self.tau[0]+1j*self.tau[1])
		self.clA()

	def clA(self):
		self.__A1_2body = 0
		self.__A2_2body = 0
		self.__A12_2body= 0

	def k1_basis(self, k1, n):
		out = np.zeros(self.__n_up+self.__n_down,dtype=np.uint8)
		cfl.copy_basis(out, n)
		return out

	@property
	def geometry(self):
		gd = {}
		gd['n_up'] = self.__n_up
		gd['n_down'] = self.__n_down
		gd['n_phi'] = self.__n_phi
		gd['lx'] = self.__lx
		gd['ly'] = self.__ly
		gd['h'] = self.__lx*self.tau[1]
		gd['area'] = self.__area
		gd['tau'] = self.tau
		gd['theta'] = self.__theta
		gd['ratio'] = 1.0/self.tau[1]
		gd['level'] = self.__level
		gd['num_thread'] = self.num_thread
		return gd

	@property
	def n_up(self):
		return self.__n_up

	@property
	def n_down(self):
		return self.__n_down

	@property
	def n_phi(self):
		return self.__n_phi

	@property
	def inter(self):
		return self.__inter

	@property
	def k1_dim(self):
		return cfl.get_k1_dim()

	def setinteraction(self, intra1, intra2, inter, calQ=False):
		if not isinstance(inter, interaction.interaction):
			raise ValueError("Variable inter is not an interaction instance!")
		if not isinstance(intra1, interaction.interaction):
			raise ValueError("Variable intra1 is not an interaction instance!")
		if not isinstance(intra2, interaction.interaction):
			raise ValueError("Variable intra2 is not an interaction instance!")
		self.__inter = inter
		self.__intra1 = intra1
		self.__intra2 = intra2
		self.clA()
		if not calQ:
			return
		geom = self.geometry
		for ele in self.__intra1.twobody:
			self.__A1_2body += getattr(interaction, 'twobody'+ele['type'])(ele, geom)
		for ele in self.__intra2.twobody:
			self.__A2_2body += getattr(interaction, 'twobody'+ele['type'])(ele, geom)
		for ele in self.__inter.twobody:
			self.__A12_2body += getattr(interaction, 'twobody'+ele['type'])(ele, geom)
		return [self.__A1_2body, self.__A2_2body, self.__A12_2body]

	def copy(self):
		return copy.deepcopy(self)

	#def k1_space(self, k1):
	#	if 0 > k1 or k1 >= self.__n_phi:
	#		raise ValueError("k1 = ", k1,"N_phi = ", self.__n_phi)
	#	path = self.basepath
	#	return np.fromfile(path + str(k1) + "_" + "base" + ".bin",dtype=np.ubyte).reshape(-1,self.n)

	#def gaugetau(self, k1):
	#	basis = self.k1_space(k1)
	#	gauge = 1j*pi*(np.arange(self.n_phi))*np.arange(self.n_phi)/self.n_phi
	#	return np.exp(np.sum(gauge[basis],axis=1))

	def Hilbertspace(self, k1=0, k2=0):
		if 0 > k1 or k1 >= self.__n_phi:
			raise ValueError("k1 = ", k1,"N_phi = ", self.__n_phi)
		if 0 > k2 or k2 >= self.__n_phi//self.__q:
			raise ValueError("k2 = ", k2,"max = ", self.__n_phi//self.__q)
		if (self.n_phi!=cfl.n_phi()) or (self.n_up!=cfl.n_up()) or (self.n_down!=cfl.n_down()):
			print("init!")
			print((self.n_phi,cfl.n_phi()),(self.n_up,cfl.n_up()),(self.n_down,cfl.n_down()))
			cfl.init(self.n_phi, self.n_up, self.n_down, self.num_thread)
		cfl.create_k1_space(k1)
		cfl.create_k2_space(k2)
		self.__k1, self.__k2 = k1, k2

	def twobodyH(self, k1=0, k2=0, HermitianQ=0):
		if isinstance(self.__A1_2body, int) or isinstance(self.__A2_2body, int) or isinstance(self.__A12_2body, int):
			self.setinteraction(self.__intra1, self.__intra2, self.__inter, calQ=True)
		self.Hilbertspace(k1, k2)
		#cfl.setHermitian(HermitianQ)
		indptr = np.zeros(cfl.get_k2_dim()+1,dtype=indtype, order = 'C')
		nnz = cfl.get_k2_n_data(indptr)
		gc.collect()
		indices = np.zeros(nnz, dtype=indtype, order = 'C')
		data = np.zeros(nnz, dtype=np.complex128, order = 'C')
		cfl.hamiltonian_generate(self.__A1_2body, self.__A2_2body, self.__A12_2body, indptr, indices, data)
		return csr_fqhe(indptr, indices, data, HermitianQ, self, k1, k2, "twobody", self.__inter)

	def tok1(self, vec):
		if isinstance(vec, eigenvector_k1):
			raise ValueError()
		self.Hilbertspace(k1=vec.k1, k2=vec.k2)
		k1, k2 = vec.k1, vec.k2
		dim = self.k1_dim
		out = np.zeros(dim, dtype = np.complex128, order='C')
		cfl.k2_to_k1_transformation(vec.vec, out)
		gc.collect()
		return eigenvector_k1(out, vec.eig, self, k1, k2)

	def density(self, vec):
		if isinstance(vec, eigenvector_k1):
			raise ValueError()
		dim = self.k1_dim
		vec_k1 = self.tok1(vec)
		out = np.zeros(2*self.__n_phi, dtype = np.float64, order='C')
		cfl.occupation_number(vec_k1.vec, out)
		return out

	def to_3layers(self, vec):
		if isinstance(vec, eigenvector_k1):
			vec_k1 = vec
		elif isinstance(vec, eigvec):
			vec_k1 = vec.tok1
		else:
			raise ValueError()
		out = np.zeros(cfl.binomial(3*self.__n_phi, self.n), dtype = np.complex128, order='C')
		k1_dim = cfl.to_3layers_trans(vec_k1.vec, out)
		return eigenvector_k1(out[0:k1_dim], vec.eig, self, vec.k1, vec.k2)
'''

class QHFMFermionDL(FermionDL):
	def __init__(self, n, n_o, dtype=np.float64, projectionQ=0, num_thread=8):
		super().__init__(n, n_o, num_thread)
		self.__level = 0

		self.__A_2body_inter = 0
		self.__A_2body_intra = 0
		self.__cdw = np.zeros(2, dtype=np.float64)
		self.__t = 0.0
		#self.__A_3body = 0
		cfl.init_FQHE(self.n_o, self.n)
		self.__dtype = dtype
		self.__qn = {}
		self.__projectionQ = projectionQ
		cfl.set_Projection(self.__projectionQ)

	def __repr__(self):
		return str(self.geometry)

	def set_Projection(self, projectionQ):
		self.__projectionQ = projectionQ
		cfl.set_Projection(self.__projectionQ)

	def clA(self):
		self.__A_2body_inter = 0
		self.__A_2body_intra = 0
		#self.__A_3body = 0
		self.__cdw = np.zeros(2, dtype=np.float64)
		self.__t = 0.0

	@property
	def geometry(self):
		gd = {}
		gd['n'] = self.n
		gd['n_phi'] = (self.n_o-1)/2
		gd['n_o'] = self.n_o
		gd['level'] = self.__level
		gd['num_thread'] = self.num_thread
		gd['qn'] = self.__qn
		return gd

	@property
	def inter(self):
		return self.__inter

	@property
	def space_dim(self):
		return cfl.get_Lz_dim()

	def setinteraction(self, intra, inter, cdw=np.zeros(2, dtype=np.float64), t=0.0, calQ=False):
		if not isinstance(inter, interaction.interaction):
			raise ValueError("Variable inter is not an interaction instance!")
		self.__inter = inter
		self.__intra = intra
		self.__twobody_inter = inter.twobody
		self.__twobody_intra = intra.twobody
		self.__A_2body_inter = 0
		self.__A_2body_intra = 0
		self.__cdw = cdw
		self.__t = t
		if not calQ:
			return
		geom = self.geometry
		#for ele in self.__twobody:
			#self.__A_2body += getattr(interaction, 'twobody'+ele['type'])(ele, geom)
		self.__A_2body_inter = np.zeros((self.n_o,self.n_o,self.n_o,self.n_o), dtype=np.float64)
		self.__A_2body_intra = np.zeros((self.n_o,self.n_o,self.n_o,self.n_o), dtype=np.float64)
		cfl.createA2BdoywithNo(self.__twobody_inter, self.__A_2body_inter, self.n_o, self.num_thread)
		cfl.createA2BdoywithNo(self.__twobody_intra, self.__A_2body_intra, self.n_o, self.num_thread)
		return self.__A_2body_inter

	def debug(self, job="A_2body_inter"):
		print("n_o: ", cfl.B_n_o(), "\tn_e: ", cfl.B_n_e(), "\tLz: ", cfl.B_n_e() )
		return self.__A_2body_inter

	def copy(self):
		return copy.deepcopy(self)

	'''

		def Lz_space(self, Lz=0, copy=False):
			#if 0 > Lz :
			#	raise ValueError("k1 = ", k1,"N_phi = ", self.n_phi)
			cfl.create_Lz_space( round(Lz+self.n*(self.n_o-1)//2) )
			self.__Lz = Lz
			self.__Z2 = None
			if copy:
				basis = np.zeros((cfl.get_Lz_dim(), self.n), dtype=np.int8)
				cfl.copy_basis(basis)
				return basis

		def Lz_Z2_space(self, Lz=0, Z2=-1, copy=False):
			cfl.create_Lz_Z2_space( round(Lz+self.n*(self.n_o-1)//2), Z2 )
			self.__Lz = Lz
			self.__Z2 = Z2
	'''
	def Hilbertspace(self, qn={"Lz":0}):
		#if 0 > k1 or k1 >= self.n_phi:
		#	raise ValueError("k1 = ", k1,"N_phi = ", self.n_phi)
		Lz = qn.get("Lz")
		Z2 = qn.get("Z_2")
		PH = qn.get("PH")
		if (self.n_o!=cfl.n_o()) or (self.n!=cfl.n_e()):
			print("init!")
			print((self.n_o,cfl.n_o()),(self.n,cfl.n_e()))
			cfl.init_FQHE(self.n_o, self.n)
		if Z2==None and PH==None and isinstance(Lz, int):
			cfl.create_Lz_space( round(Lz+self.n*(self.n_o-1)//2) )
		elif Z2 in [1, -1] and isinstance(Lz, int) and PH==None:
			cfl.create_Lz_Z2_space( round(Lz+self.n*(self.n_o-1)//2), Z2 )
		elif PH in [1, -1] and Z2 in [1, -1] and isinstance(Lz, int):
			cfl.create_Lz_Z2PH_space( round(Lz+self.n*(self.n_o-1)//2), Z2, PH )
		else:
			raise ValueError("Quantum number is incorrect! qn = ", qn)
		self.__qn = qn.copy()


	def twobodyH(self, qn={"Lz":0}, HermitianQ=0):
		if isinstance(self.__A_2body_inter, int) or isinstance(self.__A_2body_intra, int):
			self.setinteraction(self.__intra, self.__inter, self.__cdw, self.__t, calQ=True)
		self.Hilbertspace(qn)

		Lz = qn.get("Lz")
		Z2 = qn.get("Z_2")
		PH = qn.get("PH")
		#cfl.setHermitian(HermitianQ)
		indptr = np.zeros(cfl.get_Lz_dim()+1,dtype=indtype, order = 'C')
		if Z2==None and PH==None:
			nnz = cfl.twoBodyCount(indptr, self.num_thread)
		elif PH==None:
			nnz = cfl.twoBodyZ2Count(indptr, self.num_thread)
		else:
			nnz = cfl.twoBodyZ2PHCount(indptr, self.num_thread)
			#print("hahah", indptr)
		if nnz <= indptr[-1]:
			nnz = indptr[-1]
		#print("twobodyH", nnz, "cfl.get_Lz_dim()", cfl.get_Lz_dim())
		#input()
		gc.collect()
		indices = np.zeros(nnz, dtype=indtype, order = 'C')
		data = np.zeros(nnz, dtype=self.__dtype, order = 'C')
		#print(self.__A_2body_intra)
		if Z2==None and PH==None:
			cfl.twoBodyGenerate(indptr, indices, data, self.__A_2body_intra, self.__A_2body_inter, self.__cdw, self.__t, self.num_thread)
		elif PH==None:
			cfl.twoBodyZ2Generate(indptr, indices, data, self.__A_2body_intra, self.__A_2body_inter, self.__t, self.num_thread)
		else:
			cfl.twoBodyZ2PHGenerate(indptr, indices, data, self.__A_2body_intra, self.__A_2body_inter, self.__t, 0, self.num_thread)
		return csr_fqhe(indptr, indices, data, HermitianQ, self, self.__qn, "twobody", self.__inter)

	def setDefectLine(self, defectPara):
		self.__defectPara = defectPara

	def twobodyDefectLineH(self, qn={"Lz":0}, HermitianQ=0):
		if isinstance(self.__A_2body_inter, int) or isinstance(self.__A_2body_intra, int):
			self.setinteraction(self.__intra, self.__inter, self.__cdw, self.__t, calQ=True)
		self.Hilbertspace(qn)

		defectH, defectT = self.__defectPara
		pl = np.zeros(self.n_o, dtype=np.float64)
		for l in range(self.n_o):
			pl[l] = np.real(2*np.pi*np.sum( sph_harm(0, l, 0,defectT )*defectH ))

		defectlist = np.zeros( self.n_o, dtype=np.float64)
		temp_list = np.zeros( (self.n_o, self.n_o), dtype=np.float64)
		for l in range(self.n_o):
			cfl.create_nlm_AwithNo(temp_list, l, 0, self.n_o, self.num_thread)
			defectlist += pl[l]*np.diag(temp_list)
		#print(defectlist)

		Lz = qn.get("Lz")
		Z2 = qn.get("Z_2")
		PH = qn.get("PH")
		#cfl.setHermitian(HermitianQ)
		indptr = np.zeros(cfl.get_Lz_dim()+1,dtype=indtype, order = 'C')
		if Z2==None and PH==None:
			nnz = cfl.twoBodyCount(indptr, self.num_thread)
		else:
			raise ValueError( "Unsupportted QN:", qn );
			#print("hahah", indptr)
		if nnz <= indptr[-1]:
			nnz = indptr[-1]
		#print("twobodyH", nnz, "cfl.get_Lz_dim()", cfl.get_Lz_dim())
		#input()
		gc.collect()
		indices = np.zeros(nnz, dtype=indtype, order = 'C')
		data = np.zeros(nnz, dtype=self.__dtype, order = 'C')
		#print(self.__A_2body_intra)

		cfl.defectLineGenerate(indptr, indices, data, self.__A_2body_intra, self.__A_2body_inter, self.__cdw, defectlist, self.__t, self.num_thread)

		Hcsr = csr_fqhe(indptr, indices, data, HermitianQ, self, self.__qn, "twobody", self.__inter)
		return Hcsr#csr_H_L2(Hcsr.holder(), self.L2H(qn, HermitianQ).holder(), self, qn, 0, 0)

	def L2H(self, qn={"Lz":0}, HermitianQ=0):
		self.Hilbertspace(qn)

		Lz = qn.get("Lz")
		Z2 = qn.get("Z_2")
		PH = qn.get("PH")
		#cfl.setHermitian(HermitianQ)
		indptr = np.zeros(cfl.get_Lz_dim()+1,dtype=indtype, order = 'C')
		if Z2==None and PH==None:
			nnz = cfl.L2Count(indptr, self.num_thread)
		elif PH==None:
			nnz = cfl.L2Z2Count(indptr, self.num_thread)
		else:
			nnz = cfl.L2Z2PHCount(indptr, self.num_thread)
		if nnz <= indptr[-1]:
			nnz = indptr[-1]
		gc.collect()
		indices = np.zeros(nnz, dtype=indtype, order = 'C')
		data = np.zeros(nnz, dtype=self.__dtype, order = 'C')
		if Z2==None and PH==None:
			cfl.L2Generate(indptr, indices, data, self.num_thread)
		elif PH==None:
			cfl.L2Z2Generate(indptr, indices, data, self.num_thread)
		else:
			cfl.L2Z2PHGenerate(indptr, indices, data, self.num_thread)

		return csr_fqhe(indptr, indices, data, HermitianQ, self, self.__qn, "L2", interaction.interaction())

	def setL2ab(self, a=1.0, b=0.0):
		self.__a = a
		self.__b = b

	def H2L2H(self, qn={"Lz":0}, HermitianQ=0):
		return csr_H_L2(self.twobodyH(qn, HermitianQ).holder(), self.L2H(qn, HermitianQ).holder(), self, qn, self.__a, self.__b)

	def n00AH(self, A=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}, HermitianQ=0):
		if isinstance(self.__A_2body_inter, int) or isinstance(self.__A_2body_intra, int):
			self.setinteraction(self.__intra, self.__inter, self.__cdw, self.__t, calQ=True)

		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		#cfl.setHermitian(HermitianQ)
		indptr = np.zeros(cfl.get_Lz_dim()+1,dtype=indtype, order = 'C')
		nnz = cfl.n00ACount(indptr, A, self.num_thread)
		if nnz <= indptr[-1]:
			nnz = indptr[-1]
		#print("twobodyH", nnz, "cfl.get_Lz_dim()", cfl.get_Lz_dim())
		#input()
		gc.collect()
		indices = np.zeros(nnz, dtype=indtype, order = 'C')
		data = np.zeros(nnz, dtype=self.__dtype, order = 'C')
		#print(self.__A_2body_intra)
		cfl.n00AGenerate(A, indptr, indices, data, self.num_thread)
		return csr_fqhe(indptr, indices, data, HermitianQ, self, self.__qn, "twobody", self.__inter)

	def D00ABH(self, olist=[1.0], A=np.array([[1.0,0],[0,-1]]), B=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}, HermitianQ=0):
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		#cfl.setHermitian(HermitianQ)
		indptr = np.zeros(cfl.get_Lz_dim()+1,dtype=indtype, order = 'C')
		nnz = cfl.D00ACount(indptr, A, B, self.num_thread)
		if nnz <= indptr[-1]:
			nnz = indptr[-1]
		#print("twobodyH", nnz, "cfl.get_Lz_dim()", cfl.get_Lz_dim())
		#input()
		gc.collect()
		indices = np.zeros(nnz, dtype=indtype, order = 'C')
		data = np.zeros(nnz, dtype=self.__dtype, order = 'C')
		#print(self.__A_2body_intra)

		olist_ = np.zeros(self.n_o, dtype=np.float64)
		for i in range(min(len(olist_), len(olist) )):
			olist_[i] = olist[i]
		A_list = np.zeros( (self.n_o, self.n_o, (2*self.n_o-1)), dtype=np.float64)
		cfl.createO00AwithNo(olist_, A_list, self.n_o, self.num_thread)
		#return A_list;
		cfl.D00AGenerate(A, B, A_list, indptr, indices, data, self.num_thread)
		return csr_fqhe(indptr, indices, data, HermitianQ, self, self.__qn, "twobody", self.__inter)

	def D00ABOp(self, vec, L=0, olist=[1.0], A=np.array([[1.0,0],[0,-1]]), B=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}, useL0Q=False):
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		if L<0 or L>self.n-1:
			raise ValueError("Invalid input L=", L, "\t, self.n=", self.n)

		olist_ = np.zeros(self.n_o, dtype=np.float64)
		for i in range(min(len(olist_), len(olist) )):
			olist_[i] = olist[i]
		A_list = np.zeros( (self.n_o, self.n_o, (2*self.n_o-1)), dtype=np.float64)
		if L==0 and useL0Q==False:
			cfl.createO00AwithNo(olist_, A_list, self.n_o, self.num_thread)
		else:
			cfl.createOL0AwithNo(olist_, A_list, L, self.n_o, self.num_thread)

		resu = np.zeros(vec.shape[0], dtype=vec.dtype)
		cfl.D00AOp(A, B, A_list, vec, resu, self.num_thread)
		return resu

	def nnl0ABOp(self, vec, L=0, A=np.array([[1.0,0],[0,-1]]), B=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}):
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		if L<0 or L>self.n-1:
			raise ValueError("Invalid input L=", L, "\t, self.n=", self.n)

		A_list = np.zeros( (self.n_o, self.n_o, (2*self.n_o-1)), dtype=np.float64)
		cfl.creatennL0AwithNo(A_list, L, self.n_o, self.num_thread)

		resu = np.zeros(vec.shape[0], dtype=vec.dtype)
		cfl.D00AOp(A, B, A_list, vec, resu, self.num_thread)
		return resu

	def nlm_nlm_ABOp(self, vec, L=0, olist=np.array([1]), A=np.array([[1.0,0],[0,-1]]), B=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}):
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		if L<0 or L>self.n-1:
			raise ValueError("Invalid input L=", L, "\t, self.n=", self.n)

		A_list = np.zeros( (self.n_o, self.n_o, (2*self.n_o-1)), dtype=np.float64)
		cfl.create_nlm_nlm_AwithNo(A_list, olist, L, self.n_o, self.num_thread)

		resu = np.zeros(vec.shape[0], dtype=vec.dtype)
		cfl.D00AOp(A, B, A_list, vec, resu, self.num_thread)
		return resu

	def nl1m_nl2m_ABOp(self, vec, l1=0, l2=0, m=0,  A=np.array([[1.0,0],[0,-1]]), B=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}):
		'''
			n_{l1,-m} n_{l2, m}
		'''
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		if l1<0 or l1>=self.n_o or l2<0 or l2>=self.n_o:
			raise ValueError("Invalid input l_1=", l1, "\t, l_2=", l2, "\t, self.n=", self.n)
		if abs(m)>l1 or abs(m)>l2:
			raise ValueError("Invalid input l_1=", l1, "\t, l_2=", l2, "\t, m=", m)

		A_list = np.zeros( (self.n_o, self.n_o, (2*self.n_o-1)), dtype=np.float64)
		cfl.create_nl1m_nl2m_AwithNo(A_list, l1, l2, m, self.n_o, self.num_thread)
		#return A_list

		resu = np.zeros(vec.shape[0], dtype=vec.dtype)
		cfl.D00AOp(A, B, A_list, vec, resu, self.num_thread)
		return resu

	def nl0AOp(self, vec, l=0, A=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}):
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})

		resu = np.zeros(vec.shape[0], dtype=vec.dtype)
		cfl.nl0AOp(A, vec, resu, l, self.num_thread)
		return resu

	def nlmAOp(self, vec, l=0, m=0, A=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}):
		if l<0 or l>=self.n_o:
			raise ValueError("Invalid angular momentum (l, No) = ", (l, self.n_o) )
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})

		cfl.init_FQHE_B(self.n_o, self.n)
		cfl.create_Lz_space_B( round(Lz+m+self.n*(self.n_o-1)//2) )
		#print(n_new, cfl.get_B_Lz_dim())

		A_list = np.zeros( (self.n_o, self.n_o), dtype=np.float64)
		cfl.create_nlm_AwithNo(A_list, l, m, self.n_o, self.num_thread)

		resu = np.zeros( cfl.get_B_Lz_dim(), dtype=vec.dtype )
		cfl.twoOpGeneral(A, A_list, vec, resu, m, self.num_thread)

		cfl.exchangeAB()
		self.Hilbertspace({"Lz":Lz+m})
		#print(n_new, cfl.get_B_Lz_dim())
		return resu

	def PairingAlmOp(self, vec, l=0, m=0, A=np.array([[1.0,0],[0,-1]]), daggerQ=0, qn={"Lz":0}):
		'''
			daggerQ:
				0: c c|>
				1: c^dagger c^dagger|>
		'''
		if l>abs( self.n_o-1 ) or l<0:
			raise ValueError("Invalid angular momentum (l, N_o) = ", (l, self.n_o) )
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})
		sign = (1 if daggerQ!=0 else -1)

		#print("self.n", self.n)
		n_new = self.n + 2*sign
		cfl.init_FQHE_B(self.n_o, n_new)
		cfl.create_Lz_space_B( round(Lz+sign*m+n_new*(self.n_o-1)//2) )
		#print(n_new, cfl.get_B_Lz_dim())

		A_list = np.zeros( (self.n_o, self.n_o), dtype=np.float64)
		cfl.create_delta_AwithNo(A_list, l, m, self.n_o, self.num_thread)

		resu = np.zeros( cfl.get_B_Lz_dim(), dtype=vec.dtype )
		cfl.pairingOp(A_list, A, vec, resu, l, m, daggerQ, self.num_thread)

		cfl.exchangeAB()
		self.reset_n(n_new)
		self.Hilbertspace({"Lz":Lz+sign*m})
		#print(n_new, cfl.get_B_Lz_dim())
		return resu

	def PairingABlmlalbOp(self, vec, l=0, m=0, la=0, lb=0, A=np.array([[1.0,0],[0,-1]]), B=np.array([[1.0,0],[0,-1]]), qn={"Lz":0}):
		if l<abs(la-lb) or l>abs(la+lb) or l<0 or la<0 or lb<0:
			raise ValueError("Invalid angular momentum (l, la, lb) = ", (l, la, lb) )
		Lz = qn.get("Lz")
		self.Hilbertspace({"Lz":Lz})

		if m!= 0:
			cfl.init_FQHE_B(self.n_o, self.n)
			cfl.create_Lz_space_B( round(Lz+m+self.n*(self.n_o-1)//2) )
			resu_dim = cfl.get_B_Lz_dim()
		else:
			resu_dim = cfl.get_Lz_dim()

		A_list = np.zeros( (self.n_o, self.n_o, self.n_o, self.n_o), dtype=np.float64)
		cfl.create_deltadelta_ABwithNo(A_list, l, m, la, lb, self.n_o, self.num_thread)

		resu = np.zeros( resu_dim, dtype=vec.dtype )
		cfl.pairingABlmOp(A, B, A_list, vec, resu, l, m, la, lb, self.num_thread)

		if m!= 0:
			cfl.exchangeAB()
			self.Hilbertspace({"Lz":Lz+m})
			#print(n_new, cfl.get_B_Lz_dim())
		return resu

	def density(self, vec):
		Lz = vec.qn.get("Lz")
		Z2 = vec.qn.get("Z_2")
		if Z2!=None:
			raise ValueError("Unsupported method for Z_2 quantum number")
		self.Hilbertspace(qn=vec.qn)
		out = np.zeros(2*self.n_o, dtype = np.float64, order='C')
		st = cfl.meausre_density(vec.vec, out)
		gc.collect()
		if st != 0:
			raise ValueError("cfl.occupation_number, st = ", st)
		return out

	def sigma_x(self, vec):
		Lz = vec.qn.get("Lz")
		Z2 = vec.qn.get("Z_2")
		if Z2!=None:
			raise ValueError("Unsupported method for Z_2 quantum number")
		self.Hilbertspace(qn=vec.qn)
		out = np.zeros(self.n_o, dtype = vec.vec.dtype, order='C')
		st = cfl.meausre_sigma_x(vec.vec, out)
		gc.collect()
		if st != 0:
			raise ValueError("cfl.occupation_number, st = ", st)
		return out

	def Z2(self, vec):
		if isinstance(vec, np.ndarray):
			out = np.zeros(1, dtype = vec.dtype, order='C')
			st = cfl.meausre_Z2(vec, out)
			if st != 0:
				raise ValueError("cfl.meausre_Z2, st = ", st)
			return out[0]
		Lz = vec.qn.get("Lz")
		Z2 = vec.qn.get("Z_2")
		if Z2!=None:
			return Z2
		self.Hilbertspace(qn=vec.qn)
		out = np.zeros(1, dtype = vec.vec.dtype, order='C')
		st = cfl.meausre_Z2(vec.vec, out)
		gc.collect()
		if st != 0:
			raise ValueError("cfl.meausre_Z2, st = ", st)
		return out[0]

	def PH(self, vec):
		if isinstance(vec, np.ndarray):
			out = np.zeros(1, dtype = vec.dtype, order='C')
			st = cfl.meausre_PH(vec, out)
			if st != 0:
				raise ValueError("cfl.meausre_PH, st = ", st)
			return out[0]
		Lz = vec.qn.get("Lz")
		Z2 = vec.qn.get("Z_2")
		self.Hilbertspace(qn=vec.qn)
		out = np.zeros(1, dtype = vec.vec.dtype, order='C')
		if Z2!=None:
			st = cfl.meausre_PH_inZ2space(vec.vec, out)
		else:
			st = cfl.meausre_PH(vec.vec, out)
		gc.collect()
		if st != 0:
			raise ValueError("cfl.meausre_PH, st = ", st)
		return out[0]

	def Z2PH_to_Lz(self, veclist):
		Lz = self.__qn.get("Lz")
		if (self.n_o!=cfl.B_n_o()) or (self.n!=cfl.B_n_e()):
			print("initB!")
			print((self.n_o,cfl.B_n_o()),(self.n,cfl.B_n_e()))
			cfl.init_FQHE_B(self.n_o, self.n)
		cfl.create_Lz_space_B( round(Lz+self.n*(self.n_o-1)//2) )

		if isinstance(veclist, np.ndarray):
			resu = np.zeros( cfl.get_B_Lz_dim(), dtype=veclist.dtype )
			cfl.LzZ2PH_to_Lz( veclist, resu, self.num_thread )
			cfl.init_FQHE_B(self.n_o, self.n)
			return resu
		elif isinstance(veclist, eigvec):
			resu = np.zeros( cfl.get_B_Lz_dim(), dtype=veclist.vec.dtype )
			cfl.LzZ2PH_to_Lz( veclist.vec, resu, self.num_thread )
			cfl.init_FQHE_B(self.n_o, self.n)
			self.Hilbertspace(qn={"Lz":self.__qn.get("Lz")})
			resu = eigvec( resu, veclist.eig, self, {"Lz":self.__qn.get("Lz")}, measured_qn=veclist.measured_qn )
			return resu
		elif isinstance(veclist, list):
			resu = []
			for vec in veclist:
				if isinstance(vec, eigvec):
					resut = np.zeros( cfl.get_B_Lz_dim(), dtype=vec.vec.dtype )
					cfl.LzZ2PH_to_Lz( vec.vec, resut, self.num_thread )
				else:
					resut = np.zeros( cfl.get_B_Lz_dim(), dtype=vec.dtype )
					cfl.LzZ2PH_to_Lz( vec, resut, self.num_thread )
				resu.append( resut )
			cfl.init_FQHE_B(self.n_o, self.n)
			self.Hilbertspace(qn={"Lz":self.__qn.get("Lz")})
			for i in range(len(resu)):
				resu[i] = eigvec( resu[i], veclist[i].eig, self, {"Lz":self.__qn.get("Lz")}, measured_qn=veclist[i].measured_qn  )
			return resu
		else:
			raise ValueError("Unsupported instance veclist = ", type(veclist))





class fqhe_dl_solver():
	'''
		instance(ins, FQHE)
		instance(inter, interaction)
		["twobody", "L2", "H2L2", "twobodyDefectLine"]
	'''
	def __init__(self, ins, intra1, intra2, inter, itype, qn={"Lz":0}, cdw=np.zeros(2, dtype=np.float64), t=0.0, hs=None, thetas=None):
		if not isinstance(ins, FermionDL):
			raise ValueError(ins, "is not a FermionDL isinstance!")
		self.__qn = qn.copy()
		self.__fqhe = ins
		if itype not in ["twobody", "L2", "H2L2", "twobodyDefectLine"]:
			raise ValueError(itype)
		self.__itype = itype
		if isinstance(self.__fqhe, QHFMFermionDL):
			self.__fqhe.setinteraction(intra1, inter, cdw, t, calQ=False)
		else:
			self.__fqhe.setinteraction(intra1, intra2, inter, calQ=False)
		self.__intra1 = intra1
		self.__intra2 = intra2
		self.__inter = inter
		#if hs!=None and thetas!=None:
		self.__defectPara = (hs, thetas)
		self.__fqhe.setDefectLine(self.__defectPara)

	def set_qn(self, qn={"Lz":0}):
		self.__qn = qn.copy()

	def setinteraction(self, intra1, intra2, inter, cdw=np.zeros(2, dtype=np.float64), t=0.0, calQ=False):
		if not isinstance(inter, interaction.interaction):
			raise ValueError("Variable inter is not an interaction instance!")
		self.__intra1 = intra1
		self.__intra2 = intra2
		self.__inter = inter
		cal = calQ
		if isinstance(self.__fqhe, QHFMFermionDL):
			self.__fqhe.setinteraction(intra1, inter, cdw, t, calQ=cal)
		else:
			self.__fqhe.setinteraction(intra1, intra2, inter, calQ=cal)

	@property
	def fqhe(self):
		return self.__fqhe

	@property
	def qn(self):
		return self.__qn

	@property
	def itype(self):
		return self.__itype

	@property
	def geometry(self):
		return self.__fqhe.geometry
	'''
		def eig(self, tol=10**(-15), vec=False, HermitianQ=0):
			csr = getattr(self.__fqhe, self.__itype+"H")(Lz=self.__Lz, Z2=self.__Z2, HermitianQ=HermitianQ)
			ge, gv = csr.lanczos(tol = tol)
			if vec:
				return eigvec(gv, ge, self.__fqhe, self.k1, self.k2)
			else:
				return ge
	'''
	def eigs(self, tol=10**(-15), HermitianQ=0):
		csr = getattr(self.__fqhe, self.__itype+"H")(qn=self.qn, HermitianQ=HermitianQ)
		ge = csr.lanczosIRL(tol=tol,k=4)[0]
		return ge

	def eigsv(self, k=4, tol=10**(-15), vec=True, ncv=None, HermitianQ=0, measureD=True):
		h_begin = time.time()
		csr = getattr(self.__fqhe, self.__itype+"H")(qn=self.qn, HermitianQ=HermitianQ)
		h_end = time.time()
		#print(csr, getattr(self.__fqhe, self.__itype+"H"))
		lanczos_begin = time.time()
		ge = csr.lanczosIRL(tol=tol,k=k,ncv=ncv)
		lanczos_end = time.time()

		out = []
		if vec:
			for i in range(k):
				out.append(eigvec(np.transpose(ge[1])[i], ge[0][i], self.__fqhe, self.qn))
				out[i].add_measured_qn({"energy":ge[0][i]})
			if(self.__itype=="H2L2"):
				print("eigsv: Measure: H2L2")
				re_L2 = csr.measureL2list(out)
				re_E = csr.measureHlist(out)
				#print(re_L2, re_E)
				#out = [out, re_L2, re_E]
				for i in range(len(out)):
					if measureD:
						out[i].PH
						out[i].Z2
					out[i].add_measured_qn({"energy":re_E[i],"L2":re_L2[i]})
					#print(out)
		else:
			out = ge[0]
		print("total ", time.time()-h_begin, ": Hamiltonian(", h_end-h_begin, ") lanczos(", lanczos_end-lanczos_begin, ")")
		return out

	def eigsvfull(self, k=4):
		csr = getattr(self.__fqhe, self.__itype+"H")(qn=self.qn, HermitianQ=0)
		ge = np.linalg.eigh(csr.csr.toarray())
		garg = np.argsort(ge[0])[:k]
		ge = [ge[0][garg], ge[1][:,garg]]
		out = []
		for i in range(min(k, len(ge[0]))):
			out.append(eigvec(np.transpose(ge[1])[i], ge[0][i], self.__fqhe, self.qn))
			out[i].add_measured_qn({"energy":ge[0][i]})
		if(self.__itype=="H2L2"):
			print("eigsv: Measure: H2L2")
			re_L2 = csr.measureL2list(out)
			re_E = csr.measureHlist(out)
			#out = [out, re_L2, re_E]
			for i in range(len(out)):
				out[i].PH
				out[i].Z2
				out[i].add_measured_qn({"energy":re_E[i],"L2":re_L2[i]})
		return out

	def Hamiltonian(self, HermitianQ=0):
		csr = getattr(self.__fqhe, self.__itype+"H")(qn=self.qn, HermitianQ=HermitianQ)
		return csr

	def writevec(self, veclist, path):
		qns = np.zeros((len(veclist), 5), dtype=np.float64  )
		for i in range(len(veclist)):
			qns[i,0] = veclist[i].get("L2")
			qns[i,1] = veclist[i].get("PH")
			qns[i,2] = veclist[i].get("Z_2")
			qns[i,3] = veclist[i].get("eig")
			qns[i,4] = veclist[i].get("energy")
			veclist[i].vec.tofile(path+"/vec"+str(i)+".bin")
		qns.tofile(path+"/quantum_number.bin")

	def readvec(self, path, indlist=[0]):
		qns = np.zeros((len(indlist), 5), dtype=np.float64  )
		qns = np.fromfile(path+"/quantum_number.bin").reshape(-1,5)
		out = []
		for i in indlist:
			qn_ = {"Lz":0}
			measured_qn_ = {"L2": qns[i,0], "Z_2": qns[i,2],"PH": qns[i,1], "energy": qns[i,4]}
			out.append(eigvec(np.fromfile(path+"/vec"+str(i)+".bin"), qns[i,3], self.__fqhe, qn_, measured_qn=measured_qn_))
		return out