
inline void spmv(struct my_csr A, MKL_Complex16* vec, MKL_Complex16* out, MKL_Complex16 a, MKL_Complex16 b)
{
	mkl_sparse_z_mv(A.opp, a, A.mat, A.descr, vec, b, out);
	return;
}


inline void d_spmv(struct my_csr A, double* vec, double* out, double a, double b)
{
	mkl_sparse_d_mv(A.opp, a, A.mat, A.descr, vec, b, out);
	return;
}