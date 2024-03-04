
void Lanczos_quench(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
					MKL_Complex16* vector, MKL_Complex16* out, MKL_INT dim, \
					double tol, double dt, MKL_INT HermitianQ);

int Lanczos_ground(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
					MKL_Complex16* vector, MKL_Complex16* out, MKL_INT dim, \
					double tol, MKL_INT HermitianQ);


void create_csr_mkl(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
								 MKL_INT dim, MKL_INT Hermitian);

//void spmv_mkl(MKL_Complex16* vec, MKL_Complex16* out, MKL_Complex16 a, MKL_Complex16 b)
void spmv_mkl(MKL_Complex16* vec, MKL_Complex16* out, double ar, double ai, double br, double bi);

void destroy_csr_mkl();








void create_csr_mkl_L2(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
								 MKL_INT dim, MKL_INT Hermitian);

void spmv_mkl_L2(MKL_Complex16* vec, MKL_Complex16* out, double ar, double ai, double br, double bi);

void destroy_csr_mkl_L2();






void create_csr_mkl_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
								 MKL_INT dim, MKL_INT Hermitian);

void d_spmv_mkl(double* vec, double* out, double a, double b);






void create_csr_mkl_L2_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
								 MKL_INT dim, MKL_INT Hermitian);

void d_spmv_mkl_L2(double* vec, double* out, double a, double b);