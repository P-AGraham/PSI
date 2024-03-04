
struct Lanczos create_lanczos(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
							MKL_INT dim, double tol, MKL_Complex16* initv, MKL_INT Hermitian)
							// 0 General Matrix 1 Hermitian matrix
{
	struct Lanczos lan;
	lan.A = create_sp_csr(indptr, indices, data, dim, Hermitian);
	lan.A_norm = cblas_dznrm2(indptr[dim], data, 1);
	lan.tol = tol;
	lan.dim = dim;
	lan.k = 0;
	lan.alpha = (double*)mkl_calloc((size_t) (lan.dim+2), sizeof(double), 64);
	lan.beta = (double*)mkl_calloc((size_t) (lan.dim+2), sizeof(double), 64);
	lan.beta_t = (double*)mkl_calloc((size_t) (lan.dim+2), sizeof(double), 64);
	lan.n = 0;
	lan.qn = 0;
	lan.q = (struct zvec*)mkl_calloc((size_t) (lan.dim+2), sizeof(struct zvec), 64);
	for(size_t i=0;i<lan.dim+2;++i)
	{
		lan.q[i].vec = NULL;
		lan.q[i].len = dim;
	}
	lan.qn = (200>(lan.dim+2))?(lan.dim+2):200;
	lan.q[0].vec = (MKL_Complex16*)mkl_calloc((size_t) (lan.dim*200), sizeof(MKL_Complex16), 64);
	for(size_t i=0;i<lan.qn;++i)
		lan.q[i].vec = lan.q[0].vec + lan.dim * i;
	lan.Ritz.len = dim;
	lan.Ritz.vec = (MKL_Complex16*)mkl_calloc((size_t) lan.dim, sizeof(MKL_Complex16), 64);
	cblas_zcopy(dim, initv, 1, lan.q[1].vec, 1);
	lan.eigs = (double*)mkl_malloc((size_t) (lan.dim+1) * sizeof(double), 64);
	lan.eigv = NULL;
	lan.quench_tol = 0;
	return lan;
}

void destroy_lanczos(struct Lanczos* lan)
{
	destroy_sp_csr(lan->A);
	lan->A_norm = 0;
	lan->tol = 1;
	lan->dim = 0;
	if(lan->alpha!=NULL)
		mkl_free(lan->alpha);
	lan->alpha = NULL;
	if(lan->beta!=NULL)
		mkl_free(lan->beta);
	lan->beta = NULL;
	if(lan->beta_t!=NULL)
		mkl_free(lan->beta_t);
	lan->beta_t = NULL;
	lan->n = 0;
	lan->k = 0;
	if(lan->q[0].vec!=NULL)
		mkl_free(lan->q[0].vec);
	for(size_t i=0;i<lan->dim+2;++i)
	{
		lan->q[i].vec = NULL;
		lan->q[i].len = 0;
	}
	if(lan->q!=NULL)
		mkl_free(lan->q);
	lan->q = NULL;
	lan->qn = 0;
	if(lan->Ritz.vec!=NULL)
		mkl_free(lan->Ritz.vec);
	lan->Ritz.vec = NULL;
	lan->Ritz.len = 0;
	if(lan->eigs!=NULL)
		mkl_free(lan->eigs);
	lan->eigs = NULL;
	if(lan->eigv!=NULL)
		mkl_free(lan->eigv);
	lan->eigv = NULL;
	lan->quench_tol = 0;
	return;
}

int Lanczos_realloc(struct Lanczos* lan)
{
	MKL_INT temp = (100+lan->qn>lan->dim+2)?(lan->dim+2):(100+lan->qn);
	MKL_Complex16* ptr = (MKL_Complex16*)mkl_realloc(lan->q[0].vec ,(size_t) (lan->dim*temp) * sizeof(MKL_Complex16));
	if(ptr==NULL)
	{	
		printf("Lanczos_realloc failed!\n");
		return -2;
	}
	lan->qn = temp;
	lan->q[0].vec = ptr;
	for(size_t i=0;i<lan->qn;++i)
		lan->q[i].vec = lan->q[0].vec + lan->dim * i;
	return 0;
}

/*
struct Lanczos
{
	struct my_csr A;
	double A_norm;			//the norm of matrix A
	double tol;				//tolerance
	MKL_INT dim;			//dimension of the original matrix
	double* alpha;			//digonal elements of tridiagonal matrix
	double* beta;			//subdigonal elements of tridiagonal matrix
	double* beta_t;
	MKL_INT n;				//dimension of the tridiagonal matrix or steps
	MKL_INT k;				//the number of converged eigenvalues
	struct zvec* q;			//Lanczos vectors & q[0]=0 q[1] is the initial vector
	struct zvec Ritz;		//Ritz vectors used to orthogonalization
	double* eigs;			//eigenvalues of tridiagonal matrix
	double* eigv;			//eigenvector of tridiagonal matrix
};
*/

#include "LanczosReO/LanczosReothogonalization.c"