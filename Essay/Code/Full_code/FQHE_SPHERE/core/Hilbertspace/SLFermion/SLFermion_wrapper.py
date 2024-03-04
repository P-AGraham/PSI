# coding: utf-8
import numpy as np
from ctypes import *
import os

cf_fqhe = CDLL(os.path.dirname(__file__) + "/fermionSL.so")

#--------------------------------------------------------------------------------------------------------------------

#if not temp.flags['C_CONTIGUOUS']: 
#    temp = np.ascontiguous(temp, dtype=temp.dtype) # 如果不是C连续的内存，必须强制转换 
#temp_ptr = cast(temp.ctypes.data, POINTER(c_int)) # a pointer to c-array
#dll.anyfunc(teamp_ptr)

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<cFQHE_fermion.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

# int init_FQHE(size_t N_o, size_t N_e);
cf_fqhe.init_FQHE.argtypes = [c_size_t, c_size_t]
init_FQHE = cf_fqhe.init_FQHE # destroy or initial the FQHE-OBJECT

#int n_o(){return (int)A.N_o;}
#int n_e(){return (int)A.N_e;}
n_o = cf_fqhe.get_n_o
n_e = cf_fqhe.get_n_e

#void set_Hermitian(MKL_INT HermitianQ)
#MKL_INT get_Hermitian(){return A.HermitianQ;}
setHermitian = cf_fqhe.set_Hermitian
getHermitian = cf_fqhe.get_Hermitian

# int create_Lz_space(size_t L_z)
#cf_fqhe.create_k1_space.argtypes = [c_size_t, ]
create_Lz_space = cf_fqhe.create_Lz_space # bytes(path + chr(0) * 20, encoding = 'utf8')

# size_t get_Lz_dim()
get_Lz_dim = cf_fqhe.get_Lz_dim # return the dimention of Lz-Hilbertspace

#int createA2Bdoy(double *pseudo, double *A_list2, int num_thread)
#cf_fqhe.createA2Bdoy.argtypes = [
#			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
#			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
#			c_int]
#createA2Bdoy = cf_fqhe.createA2Bdoy

# int createA2BdoywithNo(double *pseudo, double *A_list2, int N_o, int num_thread)
cf_fqhe.createA2BdoywithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			c_int, c_int]
createA2BdoywithNo = cf_fqhe.createA2BdoywithNo

#void copy_basis(char* basis)
cf_fqhe.copy_basis.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.int8, ndim=2,flags='C_CONTIGUOUS')]
copy_basis = cf_fqhe.copy_basis

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/Combine.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/Htwobody.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#size_t twoBodyCount(MKL_INT *indptr, int num_thread);
cf_fqhe.twoBodyCount.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1,flags='C_CONTIGUOUS')]
twoBodyCount = cf_fqhe.twoBodyCount

#size_t twoBodyGenerate(complex double *A_list, MKL_INT *indptr, \
#							MKL_INT *indices, complex double* data, int num_thread);
cf_fqhe.twoBodyGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS')]
twoBodyGenerate = cf_fqhe.twoBodyGenerate
