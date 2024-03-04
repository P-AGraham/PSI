# coding: utf-8
import numpy as np
from ctypes import *
import os

cf_fqhe = CDLL(os.path.dirname(__file__) + "/fermionDL.so")

#--------------------------------------------------------------------------------------------------------------------

#if not temp.flags['C_CONTIGUOUS']: 
#    temp = np.ascontiguous(temp, dtype=temp.dtype) # 如果不是C连续的内存，必须强制转换 
#temp_ptr = cast(temp.ctypes.data, POINTER(c_int)) # a pointer to c-array
#dll.anyfunc(teamp_ptr)
#--------------------------------------------------------------------------------------------------------------------
#---------------------------<DLFermion.cpp>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#size_t binomial(size_t m, size_t n)//C_m^n
cf_fqhe.binomial.argtypes = [c_size_t, c_size_t]
binomial = cf_fqhe.binomial

#extern "C" void init(size_t N_phi, size_t N_up, size_t N_down, int num_thread)
cf_fqhe.init.argtypes = [c_size_t, c_size_t, c_size_t, c_int]
init = cf_fqhe.init

#extern "C" void create_k1_space(size_t k1){A.create_k1_basis(k1);}
cf_fqhe.create_k1_space.argtypes = [c_size_t]
create_k1_space = cf_fqhe.create_k1_space

#extern "C" void create_k2_space(size_t k2){A.create_k2_basis(k2);}
cf_fqhe.create_k2_space.argtypes = [c_size_t]
create_k2_space = cf_fqhe.create_k2_space

#extern "C" size_t get_k1_dim(){return A.get_k1_dim();}
#extern "C" size_t get_k2_dim(){return A.get_k2_dim();}
#extern "C" size_t get_n_phi(){return A.get_n_phi();}
#extern "C" size_t get_n_up(){return A.get_n_up();}
#extern "C" size_t get_n_down(){return A.get_n_down();}
get_k1_dim = cf_fqhe.get_k1_dim
get_k2_dim = cf_fqhe.get_k2_dim
n_phi = cf_fqhe.get_n_phi
n_up = cf_fqhe.get_n_up
n_down = cf_fqhe.get_n_down

#extern "C" size_t get_k2_n_data(unsigned int* indptr)
cf_fqhe.get_k2_n_data.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1,flags='C_CONTIGUOUS')]
get_k2_n_data = cf_fqhe.get_k2_n_data

#extern "C" size_t hamiltonian_generate(complex<double>* A1_list, complex<double>* A2_list, \
#							complex<double>* A12_list, unsigned int* indptr, \
#							unsigned int* indices, complex<double>* data)
cf_fqhe.hamiltonian_generate.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1,flags='C_CONTIGUOUS')]
hamiltonian_generate = cf_fqhe.hamiltonian_generate

#extern "C" void k2_to_k1_transformation(complex<double>* ket_k2, complex<double>* ket_k1)
cf_fqhe.k2_to_k1_transformation.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS')]
k2_to_k1_transformation = cf_fqhe.k2_to_k1_transformation

#extern "C" void copy_basis(unsigned char *out, size_t i)
cf_fqhe.copy_basis.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.ubyte, ndim=1, flags='C_CONTIGUOUS'),
			c_size_t]
copy_basis = cf_fqhe.copy_basis

#extern "C" void occupation_number(complex<double>* ket, double *out)
cf_fqhe.occupation_number.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')]
occupation_number = cf_fqhe.occupation_number

#extern "C" void to_3layers_trans(complex<double>* ket_k1, complex<double> *out)
cf_fqhe.to_3layers_trans.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS')]
to_3layers_trans = cf_fqhe.to_3layers_trans