
void p(sparse_status_t st)
{
	if(st == SPARSE_STATUS_SUCCESS)
		printf("SPARSE_STATUS_SUCCESS\n");
	if(st == SPARSE_STATUS_NOT_INITIALIZED)
		printf("SPARSE_STATUS_NOT_INITIALIZED\n");
	if(st == SPARSE_STATUS_ALLOC_FAILED)
		printf("SPARSE_STATUS_ALLOC_FAILED\n");
	if(st == SPARSE_STATUS_INVALID_VALUE)
		printf("SPARSE_STATUS_INVALID_VALUE\n");
	if(st == SPARSE_STATUS_EXECUTION_FAILED)
		printf("SPARSE_STATUS_EXECUTION_FAILED\n");
	if(st == SPARSE_STATUS_INTERNAL_ERROR)
		printf("SPARSE_STATUS_INTERNAL_ERROR\n");
	if(st == SPARSE_STATUS_NOT_SUPPORTED)
		printf("SPARSE_STATUS_NOT_SUPPORTED\n");
	return;
}

struct my_csr create_sp_csr(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
								 MKL_INT dim, MKL_INT Hermitian)
{	
	struct my_csr A;
	A.opp = SPARSE_OPERATION_NON_TRANSPOSE;
	if(Hermitian==0)
	{
		A.descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	}
	else
	{
		A.descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
		A.descr.mode = SPARSE_FILL_MODE_UPPER;
		A.descr.diag = SPARSE_DIAG_NON_UNIT;
	}
	//sparse_memory_usage_t policy = (AllocQ==0)?SPARSE_MEMORY_NONE:SPARSE_MEMORY_AGGRESSIVE;
	sparse_status_t st = mkl_sparse_z_create_csr(&(A.mat), SPARSE_INDEX_BASE_ZERO, dim, dim, indptr, indptr+1, indices, data);
	printf("create:\t\t");p(st);
	//printf("%ld\t%ld\n", Hermitian,dim);
	//printf("%lu\n", sizeof(MKL_Complex16));
	//for(size_t i=0;i<dim;++i)
	//	printf("%lf\n",creal(vec[i]));
	//st = mkl_sparse_set_mv_hint(A.mat, A.opp, A.descr, 100);
	//printf("set hint:\t");p(st);
	//st = mkl_sparse_set_dotmv_hint(A.mat, A.opp, A.descr, 100);
	//printf("set hint:\t");p(st);
	//st = mkl_sparse_set_memory_hint(A.mat, policy);
	//printf("memory hint:\t");p(st);
	st = mkl_sparse_optimize(A.mat);
	printf("optimize:\t");p(st);
	return A;
}

void destroy_sp_csr(struct my_csr A)
{
	sparse_status_t st = mkl_sparse_destroy(A.mat);
	printf("destroy:\t");p(st);
	return;
}

struct my_csr create_sp_csr_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
								 MKL_INT dim, MKL_INT Hermitian)
{	
	struct my_csr A;
	A.opp = SPARSE_OPERATION_NON_TRANSPOSE;
	if(Hermitian==0)
	{
		A.descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	}
	else
	{
		A.descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
		A.descr.mode = SPARSE_FILL_MODE_UPPER;
		A.descr.diag = SPARSE_DIAG_NON_UNIT;
	}
	//sparse_memory_usage_t policy = (AllocQ==0)?SPARSE_MEMORY_NONE:SPARSE_MEMORY_AGGRESSIVE;
	sparse_status_t st = mkl_sparse_d_create_csr(&(A.mat), SPARSE_INDEX_BASE_ZERO, dim, dim, indptr, indptr+1, indices, data);
	printf("create:\t\t");p(st);
	//printf("%ld\t%ld\n", Hermitian,dim);
	//printf("%lu\n", sizeof(MKL_Complex16));
	//for(size_t i=0;i<dim;++i)
	//	printf("%lf\n",creal(vec[i]));
	//st = mkl_sparse_set_mv_hint(A.mat, A.opp, A.descr, 100);
	//printf("set hint:\t");p(st);
	//st = mkl_sparse_set_dotmv_hint(A.mat, A.opp, A.descr, 100);
	//printf("set hint:\t");p(st);
	//st = mkl_sparse_set_memory_hint(A.mat, policy);
	//printf("memory hint:\t");p(st);
	st = mkl_sparse_optimize(A.mat);
	printf("optimize:\t");p(st);
	return A;
}

#include "Level2/mv.c"