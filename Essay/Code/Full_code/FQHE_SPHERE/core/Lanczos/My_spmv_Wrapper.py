# coding: utf-8
import numpy as np
import sys
import os
from ctypes import *
Lanczos = CDLL(os.path.dirname(__file__) + "/My_spmv.so")

indtype = np.uint64

#void my_spmv_d(double* data, size_t* indices, size_t* indptr, double* vector, double* out, double a, double b, size_t dim, size_t n_threads);
Lanczos.my_spmv_d.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=indtype, ndim=1,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
            c_double, c_double, c_size_t, c_size_t]

def my_spmv(indptr, indices, data, vec, out, a, b, num_thread, dtype=np.complex128):
    if dtype==np.complex128:
        raise ValueError("Unsupportted dtype: ", dtype)
        #Lanczos.my_spmv_z(data, indices, indptr, vec, out, a, b, indptr.shape[0]-1, num_thread)
    elif dtype==np.float64:
        Lanczos.my_spmv_d(data, indices, indptr, vec, out, a, b, indptr.shape[0]-1, num_thread)

