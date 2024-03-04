# coding: utf-8
import numpy as np
import sys
import os
from ctypes import *
Lanczos = CDLL(os.path.dirname(__file__) + "/Lanczos.so")

indtype = np.uint64
inttype = c_ulong

#--------------------------------------------------------------------------------------------------------------------

#if not temp.flags['C_CONTIGUOUS']: 
#    temp = np.ascontiguous(temp, dtype=temp.dtype) # 如果不是C连续的内存，必须强制转换 
#temp_ptr = cast(temp.ctypes.data, POINTER(c_int)) # a pointer to c-array
#dll.anyfunc(teamp_ptr)

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<lanczos_quench.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#void Lanczos_quench(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
#					MKL_Complex16* vector, MKL_Complex16* out, MKL_INT dim, \
#					double tol, double dt, MKL_INT HermitianQ)
Lanczos.Lanczos_quench.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			inttype, c_double, c_double, inttype]
Lanczos_quench = Lanczos.Lanczos_quench

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<lanczos.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#int Lanczos_ground(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
#					MKL_Complex16* vector, MKL_Complex16* out, MKL_INT dim, \
#					double tol, MKL_INT HermitianQ)
Lanczos.Lanczos_ground.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			inttype, c_double, inttype]
Lanczos_ground = Lanczos.Lanczos_ground

#int Lanczos_ground_tri(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
#					MKL_Complex16* vector, MKL_Complex16* out, double *alpha_o, \
#					double *beta_o, MKL_INT dim, \
#					double tol, MKL_INT HermitianQ)
Lanczos.Lanczos_ground_tri.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			inttype, c_double, inttype]
Lanczos_ground_tri = Lanczos.Lanczos_ground_tri


#--------------------------------------------------------------------------------------------------------------------
#---------------------------<mkl_spmv.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#void create_csr_mkl(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
#								 MKL_INT dim, MKL_INT Hermitian)
Lanczos.create_csr_mkl.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			inttype, inttype]
#void create_csr_mkl_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
#								 MKL_INT dim, MKL_INT Hermitian)
Lanczos.create_csr_mkl_d.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			inttype, inttype]
def create_csr_mkl(indptr, indices, data, dim, Hermitian):
	if data.dtype==np.complex128:
		Lanczos.create_csr_mkl(indptr, indices, data, dim, Hermitian)
	elif data.dtype==np.float64:
		Lanczos.create_csr_mkl_d(indptr, indices, data, dim, Hermitian)
	return 



#void spmv_mkl(MKL_Complex16* vec, MKL_Complex16* out, double ar, double ai, double br, double bi)
Lanczos.spmv_mkl.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_double, c_double, c_double]
#void d_spmv_mkl(double* vec, MKL_Complex16* out, double a, double b)
Lanczos.d_spmv_mkl.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_double]
def spmv_mkl(vec, out, a, b, dtype=np.complex128):
	if dtype==np.complex128:
		Lanczos.spmv_mkl(vec, out, a.real, a.imag, b.real, b.imag)
	elif dtype==np.float64:
		Lanczos.d_spmv_mkl(vec, out, a.real, b.real)


#void destroy_csr_mkl()
destroy_csr_mkl = Lanczos.destroy_csr_mkl

#------------------------------------------------------------------------------------------------

#void create_csr_mkl(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
#								 MKL_INT dim, MKL_INT Hermitian)
Lanczos.create_csr_mkl_L2.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			inttype, inttype]
#void create_csr_mkl_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
#								 MKL_INT dim, MKL_INT Hermitian)
Lanczos.create_csr_mkl_L2_d.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			inttype, inttype]
def create_csr_mkl_L2(indptr, indices, data, dim, Hermitian):
	if data.dtype==np.complex128:
		Lanczos.create_csr_mkl_L2(indptr, indices, data, dim, Hermitian)
	elif data.dtype==np.float64:
		Lanczos.create_csr_mkl_L2_d(indptr, indices, data, dim, Hermitian)
	return 

#void spmv_mkl_L2(MKL_Complex16* vec, MKL_Complex16* out, double ar, double ai, double br, double bi)
Lanczos.spmv_mkl_L2.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_double, c_double, c_double]
#void d_spmv_mkl_L2(double* vec, MKL_Complex16* out, double a, double b)
Lanczos.d_spmv_mkl_L2.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_double]
def spmv_mkl_L2(vec, out, a, b, dtype=np.complex128):
	if dtype==np.complex128:
		Lanczos.spmv_mkl_L2(vec, out, a.real, a.imag, b.real, b.imag)
	elif dtype==np.float64:
		Lanczos.d_spmv_mkl_L2(vec, out, a.real, b.real)


#void destroy_csr_mkl_L2()
destroy_csr_mkl_L2 = Lanczos.destroy_csr_mkl_L2


#void my_spmv_z(MKL_Complex16* data, size_t* indices, size_t* indptr, MKL_Complex16* vector, MKL_Complex16* out, double a, double b, size_t dim, size_t n_threads);
#void my_spmv_d(double* data, size_t* indices, size_t* indptr, double* vector, double* out, double a, double b, size_t dim, size_t n_threads);
Lanczos.my_spmv_z.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_double, c_size_t, c_size_t]
Lanczos.my_spmv_d.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_double, c_size_t, c_size_t]
def my_spmv(indptr, indices, data, vec, out, a, b, num_thread, dtype=np.complex128):
	if dtype==np.complex128:
		Lanczos.my_spmv_z(data, indices, indptr, vec, out, a, b, indptr.shape[0]-1, num_thread)
	elif dtype==np.float64:
		Lanczos.my_spmv_d(data, indices, indptr, vec, out, a, b, indptr.shape[0]-1, num_thread)