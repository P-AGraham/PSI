# coding: utf-8
import numpy as np
import gc, os
from scipy import sparse

#from .Lanczos import Lanczos_Wrapper as lanczos
from .Lanczos import My_spmv_Wrapper as lanczos
from ..general import FQHE, abpath, mkdir
from .eigenvector import eigenvector as eigvec

indtype = np.uint64

class LanczosError(Exception):
		def __init__(self, value):
			self.value = value
		def __str__(self):
			return repr(self.value)

class csr_holder():
	def __init__(self, indptr, indices, data, HermitianQ):
		if (indices.dtype != indtype) and (indptr.dtype != indtype):
			raise ValueError("indices.dtype = ",indices.dtype,
						"or indptr.dtype = ",indptr.dtype,"is not np.uint32")
		if not isinstance(HermitianQ, int):
			raise ValueError("HermitianQ must be an integer!")
		if data.dtype not in [np.complex128, np.float64]:
			raise ValueError("data.dtype = ", data.dtype, "is not np.complex128 or np.float64")
		self.__data = data
		self.__indices = indices
		self.__indptr = indptr
		self.__dtype = data.dtype
		self.__shape = (indptr.shape[0]-1,indptr.shape[0]-1)
		self.__HermitianQ = HermitianQ
		if not self.__data.flags['C_CONTIGUOUS']:
			self.__data = np.ascontiguousarray(self.__data, dtype = self.__data.dtype)
		if not self.__indices.flags['C_CONTIGUOUS']:
			self.__indices = np.ascontiguousarray(self.__indices, dtype = self.__indices.dtype)
		if not self.__indptr.flags['C_CONTIGUOUS']:
			self.__indptr = np.ascontiguousarray(self.__indptr, dtype = self.__indptr.dtype)

	@property
	def csr(self):
		if self.HermitianQ != 0:
			raise ValueError("HermitianQ = ",self.HermitianQ,"is not 0.")
		return sparse.csr_matrix((self.data, self.indices, self.indptr), shape=self.shape)

	@property
	def HermitianQ(self):
		return self.__HermitianQ

	@property
	def dtype(self):
		return self.__dtype

	@property
	def shape(self):
		return self.__shape

	@property
	def indptr(self):
		return self.__indptr
	
	@property
	def indices(self):
		return self.__indices
	
	@property
	def data(self):
		return self.__data


class csr():
	def __init__(self, indptr, indices, data, HermitianQ, num_thread=8, method="my_spmv"):
		if (indices.dtype != indtype) and (indptr.dtype != indtype):
			raise ValueError("indices.dtype = ",indices.dtype,
						"or indptr.dtype = ",indptr.dtype,"is not np.uint32")
		if not isinstance(HermitianQ, int):
			raise ValueError("HermitianQ must be an integer!")
		if data.dtype not in [np.complex128, np.float64]:
			raise ValueError("data.dtype = ", data.dtype, "is not np.complex128 or np.float64")
		self.__data = data
		self.__indices = indices
		self.__indptr = indptr
		self.__dtype = data.dtype
		self.__shape = (indptr.shape[0]-1,indptr.shape[0]-1)
		self.__HermitianQ = HermitianQ
		if not self.__data.flags['C_CONTIGUOUS']:
			self.__data = np.ascontiguousarray(self.__data, dtype = self.__data.dtype)
		if not self.__indices.flags['C_CONTIGUOUS']:
			self.__indices = np.ascontiguousarray(self.__indices, dtype = self.__indices.dtype)
		if not self.__indptr.flags['C_CONTIGUOUS']:
			self.__indptr = np.ascontiguousarray(self.__indptr, dtype = self.__indptr.dtype)
		self.__mat_n = 0
		self.__method = method #// "my_spmv" or "MKL"
		self.__num_thread = num_thread

	def holder(self):
		return csr_holder(self.__indptr, self.__indices, self.__data, self.__HermitianQ)

	def create(self):
		if self.__method == "my_spmv":
			return
		lanczos.create_csr_mkl(self.indptr,self.indices,self.data,
							self.shape[0], self.HermitianQ)
		gc.collect()

	def destroy(self):
		if self.__method == "my_spmv":
			return
		lanczos.destroy_csr_mkl()
		gc.collect()

	# 请勿外部直接调用！！！ 外部调用请使用 matvecout(self, vector)
	def matvec(self, vector):
		out = np.zeros(self.shape[0], dtype=self.__dtype, order = 'C')
		if self.__method == "my_spmv":
			lanczos.my_spmv(self.__indptr, self.__indices, self.__data, vector, out, 1.0, 0.0, self.__num_thread, dtype=self.__dtype)
		else:
			lanczos.spmv_mkl(vector, out, 1.0, 0.0, dtype=self.__dtype)
		self.__mat_n += 1
		return out

	def matvecout(self, vector):
		if self.__method == "my_spmv":
			return self.matvec(vector)
		out = np.zeros(self.shape[0], dtype=self.__dtype, order = 'C')
		self.create()
		lanczos.spmv_mkl(vector, out, 1.0, 0.0, dtype=self.__dtype)
		self.destroy()
		gc.collect()
		return out

	def lanczosIRL(self, k=4, tol=10**(-14), ncv=None):
		if self.__dtype==np.complex128:
			vec = np.random.random(self.shape[0]) + 1j*np.random.random(self.shape[0])
		elif self.__dtype==np.float64:
			vec = np.random.random(self.shape[0])
		vec /= np.linalg.norm(vec)
		self.create()
		self.__mat_n = 0
		out = sparse.linalg.eigsh(
							sparse.linalg.aslinearoperator(self), which="SA", k=k, tol=tol, v0=vec, ncv=ncv)
		self.destroy()
		gc.collect()
		print("matvec(self, vector): ", self.__mat_n, " times")
		return out
	'''
		def lanczos(self, tol=10**(-15), full_diag=True, lanczosIRL=False):
			vec = np.random.random(self.shape[0]) + 1j*np.random.random(self.shape[0])
			vec /= np.linalg.norm(vec)
			out = np.zeros(self.shape[0]+1, dtype=self.__dtype, order = 'C')
			eig = lanczos.Lanczos_ground(self.indptr, self.indices, self.data, 
								vec, out, self.shape[0], tol, self.HermitianQ)
			if eig != 0:
				print("return value: ", eig, "min_eig = ", out[-1].real)
				if self.shape[0]<=100 and full_diag:
					print("Using Full-diagonalization")
					eig = np.linalg.eigh(self.csr.toarray())
				elif self.shape[0]>100 and lanczosIRL:
					print("Using IRL-lanczos-diagonalization")
					eig = self.lanczosIRL(tol=tol)
				else:
					raise LanczosError("Lanczos is not converged!")
				out[-1] = eig[0][0]
				out[0: self.shape[0]] = np.transpose(eig[1])[0]
			gc.collect()
			return out[-1].real, out[0:self.shape[0]]
	'''

	@property
	def csr(self):
		if self.HermitianQ != 0:
			raise ValueError("HermitianQ = ",self.HermitianQ,"is not 0.")
		return sparse.csr_matrix((self.data, self.indices, self.indptr), shape=self.shape)

	'''
	def quench(self, vector, dt, tol=10**(-15)):
		if (vector.shape[0]!=self.shape[0]) or (vector.dtype!=np.complex128):
			raise ValueError("vector.shape[0] = ", vector.shape[0], "is not", self.shape[0],
							"vector.dtype = ", vector.dtype, "is not np.complex128")
		out = np.zeros(self.shape[0], dtype=np.complex128, order = 'C')
		lanczos.Lanczos_quench(self.indptr, self.indices, self.data, vector, \
								out, self.shape[0], tol, dt, self.HermitianQ)
		gc.collect()
		return out
	'''

	def expectationenergy(self, vec):
		if vec.shape[0] != self.shape[0] or vec.dtype != self.__dtype:
			raise ValueError()
		self.create()
		out = self.matvec(vec)
		self.destroy()
		return np.vdot(out, vec).real

	@property
	def HermitianQ(self):
		return self.__HermitianQ

	@property
	def dtype(self):
		return self.__dtype

	@property
	def shape(self):
		return self.__shape

	@property
	def indptr(self):
		return self.__indptr
	
	@property
	def indices(self):
		return self.__indices
	
	@property
	def data(self):
		return self.__data


class csr_fqhe(csr):
	def __init__(self, indptr, indices, data, HermitianQ, ins, qn, itype, inter, method="my_spmv"):
		super().__init__(indptr, indices, data, HermitianQ, ins.num_thread, method)
		if not isinstance(ins, FQHE):
			raise ValueError(ins, "is not a FQHE isinstance!")
		self.__qn = qn.copy()
		self.__fqhe = ins.copy()
		if itype not in ["twobody", "L2"]:
			raise ValueError(itype)
		self.__itype = itype
		self.__twobody = inter.twobody.copy()
		#self.__threebody = inter.threebody.copy()

	@property
	def HermitianTest(self):
		return sparse.linalg.norm(self.csr-self.csr.getH())

	@property
	def fqhe(self):
		return self.__fqhe

	@property
	def qn(self):
		return self.__qn

	@property
	def itype(self):
		return self.__itype
	'''
	def eig(self, tol=10**(-15), file=False, vec=True):
		calQ = True
		path1 = abpath + self.filename
		path = path1 + "k1_" + str(self.k1) + "k2_" + str(self.k2)
		eig_path = path + "_eigenvalue.bin"
		if not os.path.exists(eig_path):
			calQ = False
		if vec:
			eiv_path = path + "_eigenvector.bin"
			if not os.path.exists(eiv_path):
				calQ = False
		if not calQ:
			toll = tol
			ge, gv = self.lanczos(tol = toll)
			if file:
				mkdir(path1)
				gv.tofile(path + "_eigenvector.bin")
				np.array([ge]).tofile(eig_path)
			if vec:
				return ge, gv
			else:
				return ge
		else:
			ge = np.fromfile(eig_path, dtype=np.float64)[0]
			if vec:
				gv = np.fromfile(path + "_eigenvector.bin", dtype=np.complex128)
				return ge, gv
			return ge

	def eigs(self, tol=10**(-15), file=True):
		calQ = True
		path = abpath + self.filename
		eig_path = path + "_eigenvalue_" + str(4) + ".bin"
		if not os.path.exists(eig_path):
			toll = tol
			ge = self.lanczosIRL(tol=toll,k=4)[0]
			if file:
				ge.tofile(eig_path)
		else:
			ge = np.fromfile(eig_path, dtype=np.float64)
		return ge

	def quench(self, vector, dt, tol=10**(-15)):
		out = super().quench(vector.vec, dt, tol=tol)
		return eigvec(out, np.NaN, self.__fqhe, self.k1, self.k2)

	def matmulvec(self, vector):
		out = super().matvecout(vector.vec)
		return eigvec(out, np.NaN, self.__fqhe, self.k1, self.k2)
	'''

class csr_H_L2():
	def __init__(self, csr_holder_H, csr_holder_L2, ins, qn, a, b, method="my_spmv"):
		if not isinstance(csr_holder_H, csr_holder):
			raise ValueError(csr_holder_H, "is not a csr_holder isinstance!")
		if not isinstance(csr_holder_L2, csr_holder):
			raise ValueError( csr_holder_L2, "is not a csr_holder isinstance!")
		if not isinstance(ins, FQHE):
			raise ValueError(ins, "is not a FQHE isinstance!")
		if csr_holder_H.dtype != csr_holder_L2.dtype:
			raise ValueError(csr_holder_H.dtype, " != ", csr_holder_L2.dtype)
		if csr_holder_H.shape != csr_holder_L2.shape:
			raise ValueError(csr_holder_H.shape, " != ", csr_holder_L2.shape)
		self.__H = csr_holder_H
		self.__L2 = csr_holder_L2
		self.__qn = qn.copy()
		self.__fqhe = ins.copy()
		self.__itype = "H2+a*(L2-b)^2"
		self.__a = a
		self.__b = b
		if self.__H.dtype == self.__L2.dtype:
			self.__dtype = self.__H.dtype
		else:
			raise ValueError("self.__H.dtype: ", self.__H.dtype, " != ", "self.__L2.dtype: ", self.__L2.dtype)
		self.__mat_n = 0
		self.__method = method #// "my_spmv" or "MKL"

	@property
	def dtype(self):
		return self.__H.dtype

	@property
	def shape(self):
		return self.__L2.shape

	@property
	def HermitianTest(self):
		return (sparse.linalg.norm(self.__H.csr-self.__H.csr.getH()),sparse.linalg.norm(self.__L2.csr-self.__L2.csr.getH()))

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
	def H2(self):
		return self.__H

	@property
	def L2(self):
		return self.__L2

	@property
	def csr(self):
		return self.H2.csr + self.__a * (self.L2.csr-self.__b)@(self.L2.csr-self.__b)

	def checktype(self, type):
		if type not in ["H2+a*(L2-b)^2", "H2", "L2"]:
			raise ValueError(itype)

	def create(self):#, type="H2+a*(L2-b)^2"):
		if self.__method == "my_spmv":
			return
		#self.checktype(type)
		#if type in ["H2+a*(L2-b)^2", "H2"]:
		lanczos.create_csr_mkl(self.__H.indptr,self.__H.indices,self.__H.data,
							self.__H.shape[0], self.__H.HermitianQ)
		#if type in ["H2+a*(L2-b)^2", "L2"]:
		lanczos.create_csr_mkl_L2(self.__L2.indptr,self.__L2.indices,self.__L2.data,
							self.__L2.shape[0], self.__L2.HermitianQ)
		gc.collect()

	def destroy(self):#, type="H2+a*(L2-b)^2"):
		if self.__method == "my_spmv":
			return
		#self.checktype(type)
		#if type in ["H2+a*(L2-b)^2", "H2"]:
		lanczos.destroy_csr_mkl()
		#if type in ["H2+a*(L2-b)^2", "L2"]:
		#lanczos.destroy_csr_mkl_L2()
		gc.collect()

	# 请勿外部直接调用！！！ 外部调用请使用 matvecout(self, vector)
	def matvec(self, vector):
		'''
			lanczos.spmv_mkl(vec, out, a, b)
				out = 2*H*vec + b*out
		'''
		self.__mat_n += 1
		num_thread = self.__fqhe.num_thread
		if self.__method == "my_spmv":
			if abs(self.__b)>1e-14:
				out1 = vector.copy()
				lanczos.my_spmv(self.L2.indptr, self.L2.indices, self.L2.data, vector, out1, 1.0, -self.__b, num_thread, dtype=self.__dtype)
				out = out1.copy()
				if abs(self.__a)>1e-14:
					lanczos.my_spmv(self.L2.indptr, self.L2.indices, self.L2.data, out1, out, self.__a, -self.__a*self.__b, num_thread, dtype=self.__dtype)
				lanczos.my_spmv(self.H2.indptr, self.H2.indices, self.H2.data, vector, out, 1.0, 1.0, num_thread, dtype=self.__dtype)
			else:
				out = np.zeros(self.shape[0], dtype=self.__dtype, order = 'C')
				if abs(self.__a)>1e-14:
					lanczos.my_spmv(self.L2.indptr, self.L2.indices, self.L2.data, vector, out, self.__a, 0.0, num_thread, dtype=self.__dtype)
				lanczos.my_spmv(self.H2.indptr, self.H2.indices, self.H2.data, vector, out, 1.0, 1.0, num_thread, dtype=self.__dtype)
			return out
			
		#out = np.zeros(self.shape[0], dtype=np.complex128, order = 'C')
		if abs(self.__b)>1e-14:
			out1 = vector.copy()
			lanczos.spmv_mkl_L2(vector, out1, 1.0, -self.__b, dtype=self.__dtype) # (L^2-b) |vec>
			out = out1.copy()
			if abs(self.__a)>1e-14:
				lanczos.spmv_mkl_L2(out1, out, self.__a, -self.__a*self.__b, dtype=self.__dtype) # |out> = (a*L^2-ab^)(L^2-b)|vec> = a(L^2-b)(L^2-b)|vec>
			lanczos.spmv_mkl(vector, out, 1.0, 1.0, dtype=self.__dtype) # h|vec> + |out> = h|vec> + a(L^2-b)(L^2-b)|vec>
		else:
			out = np.zeros(self.shape[0], dtype=self.__dtype, order = 'C')
			if abs(self.__a)>1e-14:
				lanczos.spmv_mkl_L2(vector, out, self.__a, 0, dtype=self.__dtype) # aL^2 |vec>
			lanczos.spmv_mkl(vector, out, 1.0, 1.0, dtype=self.__dtype) # h|vec> + |out> = h|vec> + aL^2|vec>
		return out

	def matvecout(self, vector):
		if self.__method == "my_spmv":
			return self.matvec(vector)
		out = np.zeros(self.shape[0], dtype=self.__H.dtype, order = 'C')
		self.create()
		out1 = vector.copy()
		lanczos.spmv_mkl_L2(vector, out1, 1.0, -self.__b, dtype=self.__dtype) # (L^2-b) |vec>
		out = out1.copy()
		lanczos.spmv_mkl_L2(out1, out, self.__a, -self.__a*self.__b, dtype=self.__dtype) # |out> = (a*L^2-ab^)(L^2-b)|vec> = a(L^2-b)(L^2-b)|vec>
		lanczos.spmv_mkl(vector, out, 1.0, 1.0, dtype=self.__dtype) # h|vec> + |out> = h|vec> + a(L^2-b)(L^2-b)|vec>
		return out
		self.destroy()
		gc.collect()
		return out

	def lanczosIRL(self, k=4, tol=10**(-14), ncv=None):
		if self.__dtype==np.complex128:
			vec = np.random.random(self.shape[0]) + 1j*np.random.random(self.shape[0])
		elif self.__dtype==np.float64:
			vec = np.random.random(self.shape[0])
		vec /= np.linalg.norm(vec)
		self.create()
		self.__mat_n = 0
		out = sparse.linalg.eigsh(
							sparse.linalg.aslinearoperator(self), which="SA", k=k, tol=tol, v0=vec, ncv=ncv)
		self.destroy()
		gc.collect()
		print("matvec(self, vector): ", self.__mat_n, " times")
		return out

	def measureH(self, vec):
		num_thread = self.__fqhe.num_thread
		if vec.shape[0] != self.__H.shape[0] or vec.dtype != self.__H.dtype:
			raise ValueError()
		out = np.zeros(self.shape[0], dtype=self.__H.dtype, order = 'C')
		if self.__method == "my_spmv":
			lanczos.my_spmv(self.H2.indptr, self.H2.indices, self.H2.data, vec, out, 1.0, 0.0, num_thread, dtype=self.__dtype)
		else:
			lanczos.create_csr_mkl(self.__H.indptr,self.__H.indices,self.__H.data,
								self.__H.shape[0], self.__H.HermitianQ)
			lanczos.spmv_mkl(vec, out, 1.0, 0, dtype=self.__dtype) # (L^2-b) |vec>
			lanczos.destroy_csr_mkl()
		return np.vdot(out, vec).real

	def measureHlist(self, veclist):
		num_thread = self.__fqhe.num_thread
		if self.__method == "my_spmv":
			resu = np.zeros(len(veclist), dtype=np.float64)
			for i in range(len(veclist)):
				out = np.zeros(self.shape[0], dtype=self.__H.dtype, order = 'C')
				lanczos.my_spmv(self.H2.indptr, self.H2.indices, self.H2.data, veclist[i].vec, out, 1.0, 0.0, num_thread, dtype=self.__dtype)
				resu[i] = np.vdot(out, veclist[i].vec).real
			return resu
		lanczos.create_csr_mkl(self.__H.indptr,self.__H.indices,self.__H.data,
							self.__H.shape[0], self.__H.HermitianQ)
		resu = np.zeros(len(veclist), dtype=np.float64)
		for i in range(len(veclist)):
			out = np.zeros(self.shape[0], dtype=self.__H.dtype, order = 'C')
			lanczos.spmv_mkl(veclist[i].vec, out, 1.0, 0, dtype=self.__dtype) # (L^2-b) |vec>
			resu[i] = np.vdot(out, veclist[i].vec).real
		lanczos.destroy_csr_mkl()
		return resu

	def measureL2(self, vec):
		num_thread = self.__fqhe.num_thread
		if vec.shape[0] != self.__L2[0] or vec.dtype != self.__L2.dtype:
			raise ValueError()
		out = np.zeros(self.shape[0], dtype=self.__H.dtype, order = 'C')
		if self.__method == "my_spmv":
			lanczos.my_spmv(self.L2.indptr, self.L2.indices, self.L2.data, vec, out, 1.0, 0.0, num_thread, dtype=self.__dtype)
		else:
			lanczos.create_csr_mkl_L2(self.__L2.indptr,self.__L2.indices,self.__L2.data,
								self.__L2.shape[0], self.__L2.HermitianQ)
			lanczos.spmv_mkl_L2(veclist[i].vec, out, 1.0, 0, dtype=self.__dtype) # (L^2-b) |vec>
			lanczos.destroy_csr_mkl_L2()
		return np.vdot(out, vec).real

	def measureL2list(self, veclist):
		num_thread = self.__fqhe.num_thread
		if self.__method == "my_spmv":
			resu = np.zeros(len(veclist), dtype=np.float64)
			for i in range(len(veclist)):
				out = np.zeros(self.shape[0], dtype=self.__H.dtype, order = 'C')
				lanczos.my_spmv(self.L2.indptr, self.L2.indices, self.L2.data, veclist[i].vec, out, 1.0, 0.0, num_thread, dtype=self.__dtype)
				resu[i] = np.vdot(out, veclist[i].vec).real
			return resu
		lanczos.create_csr_mkl_L2(self.__L2.indptr,self.__L2.indices,self.__L2.data,
							self.__L2.shape[0], self.__L2.HermitianQ)
		resu = np.zeros(len(veclist), dtype=np.float64)
		for i in range(len(veclist)):
			out = np.zeros(self.shape[0], dtype=self.__L2.dtype, order = 'C')
			lanczos.spmv_mkl_L2(veclist[i].vec, out, 1.0, 0, dtype=self.__dtype) # (L^2-b) |vec>
			resu[i] = np.vdot(out, veclist[i].vec).real
		lanczos.destroy_csr_mkl_L2()
		return resu

	'''
	def lanczos(self, tol=10**(-15), full_diag=True, lanczosIRL=False):
		vec = np.random.random(self.shape[0]) + 1j*np.random.random(self.shape[0])
		vec /= np.linalg.norm(vec)
		out = np.zeros(self.shape[0]+1, dtype=np.complex128, order = 'C')
		eig = lanczos.Lanczos_ground(self.indptr, self.indices, self.data, 
							vec, out, self.shape[0], tol, self.HermitianQ)
		if eig != 0:
			print("return value: ", eig, "min_eig = ", out[-1].real)
			if self.shape[0]<=100 and full_diag:
				print("Using Full-diagonalization")
				eig = np.linalg.eigh(self.csr.toarray())
			elif self.shape[0]>100 and lanczosIRL:
				print("Using IRL-lanczos-diagonalization")
				eig = self.lanczosIRL(tol=tol)
			else:
				raise LanczosError("Lanczos is not converged!")
			out[-1] = eig[0][0]
			out[0: self.shape[0]] = np.transpose(eig[1])[0]
		gc.collect()
		return out[-1].real, out[0:self.shape[0]]

	@property
	def csr(self):
		if self.HermitianQ != 0:
			raise ValueError("HermitianQ = ",self.HermitianQ,"is not 0.")
		return sparse.csr_matrix((self.data, self.indices, self.indptr), shape=self.shape)

	def quench(self, vector, dt, tol=10**(-15)):
		if (vector.shape[0]!=self.shape[0]) or (vector.dtype!=np.complex128):
			raise ValueError("vector.shape[0] = ", vector.shape[0], "is not", self.shape[0],
							"vector.dtype = ", vector.dtype, "is not np.complex128")
		out = np.zeros(self.shape[0], dtype=np.complex128, order = 'C')
		lanczos.Lanczos_quench(self.indptr, self.indices, self.data, vector, \
								out, self.shape[0], tol, dt, self.HermitianQ)
		gc.collect()
		return out

	def expectationenergy(self, vec):
		if vec.shape[0] != self.shape[0] or vec.dtype != np.complex128:
			raise ValueError()
		self.create()
		out = self.matvec(vec)
		self.destroy()
		return np.vdot(out, vec).real

	'''
