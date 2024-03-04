# coding: utf-8
import numpy as np
from ctypes import *
import os


indtype = np.uint64
gdtype = np.float64
#gdtype = np.complex128

if gdtype==np.complex128:
	cf_fqhe = CDLL(os.path.dirname(__file__) + "/QHFM_fermion_complex.so")
elif gdtype==np.float64:
	cf_fqhe = CDLL(os.path.dirname(__file__) + "/QHFM_fermion_real.so")

#--------------------------------------------------------------------------------------------------------------------

#if not temp.flags['C_CONTIGUOUS']: 
#    temp = np.ascontiguous(temp, dtype=temp.dtype) # 如果不是C连续的内存，必须强制转换 
#temp_ptr = cast(temp.ctypes.data, POINTER(c_int)) # a pointer to c-array
#dll.anyfunc(teamp_ptr)

#extern "C" double cg_py(double la, double lb, double l, double m1, double m2, double m3);
cf_fqhe.cg_py.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double]
cf_fqhe.cg_py.restype = c_double
cg_py = cf_fqhe.cg_py # destroy or initial the FQHE-OBJECT

#extern "C" int wigner_3j_py(double la, double lb, double l, double m1, double m2, double m3);
cf_fqhe.wigner_3j_py.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double]
cf_fqhe.wigner_3j_py.restype = c_double
wigner_3j_py = cf_fqhe.wigner_3j_py

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<cFQHE_fermion.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

# extern "C" int init_FQHE(size_t N_o, size_t N_e)
cf_fqhe.init_FQHE.argtypes = [c_size_t, c_size_t]
cf_fqhe.init_FQHE_B.argtypes = [c_size_t, c_size_t]
init_FQHE = cf_fqhe.init_FQHE # destroy or initial the FQHE-OBJECT
init_FQHE_B = cf_fqhe.init_FQHE_B

# extern "C" int get_n_o(){return (int)A.N_o;}
# extern "C" int get_n_e(){return (int)A.N_e;}
n_o = cf_fqhe.get_n_o
n_e = cf_fqhe.get_n_e
B_n_o = cf_fqhe.get_B_n_o
B_n_e = cf_fqhe.get_B_n_e

#void set_Hermitian(MKL_INT HermitianQ)
#MKL_INT get_Hermitian(){return A.HermitianQ;}
#setHermitian = cf_fqhe.set_Hermitian
#getHermitian = cf_fqhe.get_Hermitian

# extern "C" int create_Lz_space(size_t L_z)
#cf_fqhe.create_k1_space.argtypes = [c_size_t, ]
create_Lz_space = cf_fqhe.create_Lz_space # bytes(path + chr(0) * 20, encoding = 'utf8')
create_Lz_space_B = cf_fqhe.create_Lz_space_B


# extern "C" size_t get_Lz_dim(){return A.dim;}
cf_fqhe.get_Lz_dim.restype = c_size_t
cf_fqhe.get_B_Lz_dim.restype = c_size_t
get_Lz_dim = cf_fqhe.get_Lz_dim # return the dimention of Lz-Hilbertspace
get_B_Lz_dim = cf_fqhe.get_B_Lz_dim 

#int createA2Bdoy(double *pseudo, double *A_list2, int num_thread)
#cf_fqhe.createA2Bdoy.argtypes = [
#			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
#			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
#			c_int]
#createA2Bdoy = cf_fqhe.createA2Bdoy

# extern "C" void copy_basis(char* basis)
cf_fqhe.copy_basis.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.int8, ndim=2,flags='C_CONTIGUOUS')]
copy_basis = cf_fqhe.copy_basis

# extern "C" void set_Projection(int ProjectionQ)
cf_fqhe.set_Projection.argtypes = [c_int]
set_Projection = cf_fqhe.set_Projection

#extern "C" int createA2BdoywithNo(double *pseudo, double *A_list2, int N_o, int num_thread);
cf_fqhe.createA2BdoywithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'), c_int, c_int]
createA2BdoywithNo = cf_fqhe.createA2BdoywithNo

#extern "C" int createO00AwithNo(double *olist, double *A_list2, int N_o, int num_thread);
cf_fqhe.createO00AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=3,flags='C_CONTIGUOUS'), c_int, c_int]
createO00AwithNo = cf_fqhe.createO00AwithNo

#extern "C" int createOL0AwithNo(double *olist, double *A_list2, int L, int N_o, int num_thread);
cf_fqhe.createOL0AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=3,flags='C_CONTIGUOUS'), c_int, c_int, c_int]
createOL0AwithNo = cf_fqhe.createOL0AwithNo

#extern "C" int creatennL0AwithNo(double *A_list2, int L, int N_o, int num_thread);
cf_fqhe.creatennL0AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=3,flags='C_CONTIGUOUS'), c_int, c_int, c_int]
creatennL0AwithNo = cf_fqhe.creatennL0AwithNo

#extern "C" int create_nlm_nlm_AwithNo(double *A_list2, double *olist, int L, int N_o, int num_thread)
cf_fqhe.create_nlm_nlm_AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=3,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'), c_int, c_int, c_int]
create_nlm_nlm_AwithNo = cf_fqhe.create_nlm_nlm_AwithNo

#extern "C" int create_nl1m_nl2m_AwithNo(double *A_list2, int l_1, int l_2, int m, int N_o, int num_thread);
cf_fqhe.create_nl1m_nl2m_AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=3,flags='C_CONTIGUOUS'),
			c_int, c_int, c_int, c_int, c_int]
create_nl1m_nl2m_AwithNo = cf_fqhe.create_nl1m_nl2m_AwithNo

#extern "C" int exchangeAB()
exchangeAB = cf_fqhe.exchangeAB

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/Combine.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/Htwobody.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

# extern "C" size_t twoBodyCount(MKL_INT *indptr, int num_thread);
cf_fqhe.twoBodyCount.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.twoBodyCount.restype = c_size_t
twoBodyCount = cf_fqhe.twoBodyCount

# extern "C" size_t twoBodyGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
#					double *A1_list, double *A12_list, double *cdw, double t, int num_thread);
cf_fqhe.twoBodyGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_int]
cf_fqhe.twoBodyGenerate.restype = c_size_t
twoBodyGenerate = cf_fqhe.twoBodyGenerate


#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/twobody_Z2.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t twoBodyZ2Count(MKL_INT *indptr, int num_thread);
cf_fqhe.twoBodyZ2Count.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.twoBodyZ2Count.restype = c_size_t
twoBodyZ2Count = cf_fqhe.twoBodyZ2Count

#extern "C" size_t twoBodyZ2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
#                    double *A1_list, double *A12_list, double t, int num_thread);
cf_fqhe.twoBodyZ2Generate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			c_double, c_int]
cf_fqhe.twoBodyZ2Generate.restype = c_size_t
twoBodyZ2Generate = cf_fqhe.twoBodyZ2Generate

# extern "C" int create_Lz_Z2_space(size_t L_z, int Z_2)// Z_2: -1 or 1 or -2
cf_fqhe.create_Lz_Z2_space.argtypes = [c_size_t, c_int]
create_Lz_Z2_space = cf_fqhe.create_Lz_Z2_space


#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/twobody_Z2PH.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t twoBodyZ2PHCount(MKL_INT *indptr, int num_thread);
cf_fqhe.twoBodyZ2PHCount.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.twoBodyZ2PHCount.restype = c_size_t
twoBodyZ2PHCount = cf_fqhe.twoBodyZ2PHCount

#extern "C" size_t twoBodyZ2PHGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
#                    double *A1_list, double *A12_list, double t, int sortQ, int num_thread);
cf_fqhe.twoBodyZ2PHGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			c_double, c_int, c_int]
cf_fqhe.twoBodyZ2PHGenerate.restype = c_size_t
twoBodyZ2PHGenerate = cf_fqhe.twoBodyZ2PHGenerate

# extern "C" int create_Lz_Z2_space(size_t L_z, int Z_2, int PH)// Z_2: -1 or 1 or -2
cf_fqhe.create_Lz_Z2PH_space.argtypes = [c_size_t, c_int, c_int]
create_Lz_Z2PH_space = cf_fqhe.create_Lz_Z2PH_space

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/L2.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t L2Count(MKL_INT *indptr, int num_thread);
cf_fqhe.L2Count.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.L2Count.restype = c_size_t
L2Count = cf_fqhe.L2Count

#extern "C" size_t L2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);
cf_fqhe.L2Generate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.L2Generate.restype = c_size_t
L2Generate = cf_fqhe.L2Generate

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/L2_Z2.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t L2Z2Count(MKL_INT *indptr, int num_thread);
cf_fqhe.L2Z2Count.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.L2Z2Count.restype = c_size_t
L2Z2Count = cf_fqhe.L2Z2Count

#extern "C" size_t L2Z2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);
cf_fqhe.L2Z2Generate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.L2Z2Generate.restype = c_size_t
L2Z2Generate = cf_fqhe.L2Z2Generate

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/L2_Z2PH.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t L2Z2PHCount(MKL_INT *indptr, int num_thread);
cf_fqhe.L2Z2PHCount.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.L2Z2PHCount.restype = c_size_t
L2Z2PHCount = cf_fqhe.L2Z2PHCount

#extern "C" size_t L2Z2PHGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);
cf_fqhe.L2Z2PHGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
cf_fqhe.L2Z2PHGenerate.restype = c_size_t
L2Z2PHGenerate = cf_fqhe.L2Z2PHGenerate

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/nlm.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t n00ACount(MKL_INT *indptr, Data_dtype *A_matirx, int num_thread);
cf_fqhe.n00ACount.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'), c_int]
cf_fqhe.n00ACount.restype = c_size_t
n00ACount = cf_fqhe.n00ACount

#extern "C" size_t n00AGenerate(Data_dtype *A_matirx, MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);
cf_fqhe.n00AGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int]
cf_fqhe.n00AGenerate.restype = c_size_t
n00AGenerate = cf_fqhe.n00AGenerate


#extern "C" size_t D00ACount(MKL_INT *indptr, Data_dtype *A_matirx, Data_dtype *B_matirx, int num_thread)
cf_fqhe.D00ACount.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'), c_int]
cf_fqhe.n00ACount.restype = c_size_t
D00ACount = cf_fqhe.D00ACount

#extern "C" size_t D00AGenerate(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, 
#								MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread)
cf_fqhe.D00AGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=3,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int]
cf_fqhe.D00AGenerate.restype = c_size_t
D00AGenerate = cf_fqhe.D00AGenerate

#extern "C" void D00AOp(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, int num_thread)
cf_fqhe.D00AOp.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=3,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int]
D00AOp = cf_fqhe.D00AOp

#extern "C" void Dl0AOp(Data_dtype *A_matirx, Data_dtype *vec, Data_dtype *resu, int l, int num_thread)
cf_fqhe.nl0AOp.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int,c_int]
nl0AOp = cf_fqhe.nl0AOp

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/General_operators.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" void pairingOp(Data_dtype *A_list, Data_dtype *A_matirx, Data_dtype *vec, Data_dtype *resu, int l, int m, 
#                            int daggerQ, int num_thread)
cf_fqhe.pairingOp.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int,c_int,c_int,c_int]
pairingOp = cf_fqhe.pairingOp

#extern "C" int create_delta_AwithNo(double *A_list2, int l, int m, int N_o, int num_thread);
cf_fqhe.create_delta_AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),c_int,c_int,c_int,c_int]
create_delta_AwithNo = cf_fqhe.create_delta_AwithNo

#extern "C" void pairingABlmOp(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, 
#                                int l, int m, int la, int lb, int num_thread)
cf_fqhe.pairingABlmOp.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]
pairingABlmOp = cf_fqhe.pairingABlmOp

#extern "C" int create_deltadelta_ABwithNo(double *A_list2, int l, int m, int la, int lb, int N_o, int num_thread)
cf_fqhe.create_deltadelta_ABwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=4,flags='C_CONTIGUOUS'),
			c_int,c_int,c_int,c_int,c_int,c_int]
create_deltadelta_ABwithNo = cf_fqhe.create_deltadelta_ABwithNo

#extern "C" void twoOpGeneral(Data_dtype *A_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, int m, int num_thread)
cf_fqhe.twoOpGeneral.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),c_int,c_int]
twoOpGeneral = cf_fqhe.twoOpGeneral

#extern "C" int create_nlm_AwithNo(double *A_list2, int l, int m, int N_o, int num_thread);
cf_fqhe.create_nlm_AwithNo.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=2,flags='C_CONTIGUOUS'),
			c_int,c_int,c_int,c_int]
create_nlm_AwithNo = cf_fqhe.create_nlm_AwithNo

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Basis/density.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

# extern "C" int meausre_density(Data_dtype *ket, double *out);
cf_fqhe.meausre_density.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS')]
meausre_density = cf_fqhe.meausre_density

# extern "C" int meausre_sigma_x(Data_dtype *ket, Data_dtype *out);
cf_fqhe.meausre_sigma_x.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
meausre_sigma_x = cf_fqhe.meausre_sigma_x

# extern "C" int meausre_Z2(Data_dtype *ket, Data_dtype *out);
cf_fqhe.meausre_Z2.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
meausre_Z2 = cf_fqhe.meausre_Z2

# extern "C" int meausre_PH(Data_dtype *ket, Data_dtype *out);
cf_fqhe.meausre_PH.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
meausre_PH = cf_fqhe.meausre_PH

# extern "C" int meausre_PH_inZ2space(Data_dtype *ket, Data_dtype *out);
cf_fqhe.meausre_PH_inZ2space.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS')]
meausre_PH_inZ2space = cf_fqhe.meausre_PH_inZ2space


#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Basis/transformation.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

# extern "C" void LzZ2PH_to_Lz(const Data_dtype* LzZ2PH_ket, Data_dtype* Lz_ket, int num_thread)
cf_fqhe.LzZ2PH_to_Lz.argtypes = [
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'), c_int]
LzZ2PH_to_Lz = cf_fqhe.LzZ2PH_to_Lz

#--------------------------------------------------------------------------------------------------------------------
#---------------------------<Hamiltonian/lineDefect.c>--------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#extern "C" size_t defectLineGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, double *A1_list, double *A12_list, double #*cdw, double* defect_list, double t, int num_thread);
cf_fqhe.defectLineGenerate.argtypes = [
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=gdtype, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=4,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
			c_double, c_int]
cf_fqhe.defectLineGenerate.restype = c_size_t
defectLineGenerate = cf_fqhe.defectLineGenerate