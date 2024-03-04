
struct my_csr
{
	sparse_matrix_t mat;
	sparse_operation_t opp;
	struct matrix_descr descr;
};

struct my_csr H;
struct my_csr L2;

void p(sparse_status_t st);

struct my_csr create_sp_csr(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
								 MKL_INT dim, MKL_INT Hermitian);

struct my_csr create_sp_csr_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
								 MKL_INT dim, MKL_INT Hermitian);

void destroy_sp_csr(struct my_csr A);

//-------------------------------Level2/mv.c-----------------------------------------------------------

void spmv(struct my_csr A, MKL_Complex16* vec, MKL_Complex16* out, MKL_Complex16 a, MKL_Complex16 b);

void d_spmv(struct my_csr A, double* vec, double* out, double a, double b);