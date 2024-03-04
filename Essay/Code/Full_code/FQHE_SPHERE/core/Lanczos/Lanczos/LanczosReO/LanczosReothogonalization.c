
int LanczosIteration(struct Lanczos* lan)
{
	if(lan->k<0)
		return -1;
	
	++(lan->n);
	if(lan->q[lan->n+1].vec==NULL)
	{
		if(Lanczos_realloc(lan)==-2)
		{
			lan->k = -2;return -2;
		}
	}

	//lan->q[lan->n+1].vec = (MKL_Complex16*)mkl_malloc((size_t) lan->dim * sizeof(MKL_Complex16), 64);
	spmv(lan->A, lan->q[lan->n].vec, lan->q[lan->n+1].vec, 1, 0);
	lan->alpha[lan->n] = z_orthogonlization(&(lan->q[lan->n+1]), &(lan->q[lan->n]));//z=z-<qj,z>qj
	MKL_Complex16 a = -lan->beta[lan->n-1]; // beta_{n-1}
	cblas_zaxpy(lan->dim, &a, lan->q[lan->n-1].vec, 1, lan->q[lan->n+1].vec, 1); //z=z-beta_{n-1}*q_{j-1}

	//diagonalize the tridigonal matrix
	if(lan->eigv!=NULL)
		mkl_free(lan->eigv);
	lan->eigv = (double*)mkl_malloc((size_t) lan->n * lan->n * sizeof(double), 64);
	cblas_dcopy(lan->n+1, lan->alpha, 1, lan->eigs, 1);
	cblas_dcopy(lan->n+1, lan->beta, 1, lan->beta_t, 1);
	LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', lan->n, lan->eigs+1, lan->beta_t+1, lan->eigv, lan->n);

	//cal beta and tolerance
	lan->beta[lan->n] = z_norm(&(lan->q[lan->n+1]));
	for(size_t i=0;i<lan->n;++i)
		z_orthogonlization(&(lan->q[lan->n+1]), &(lan->q[i+1]));
	lan->beta[lan->n] = z_normalization(&(lan->q[lan->n+1]));

	//cal the quench tolerance
	if(lan->n==1)
		lan->quench_tol = 1;
	else
		lan->quench_tol *= ((double)(lan->n-1)/lan->beta[lan->n-1]);
	//printf("%d\t%e\n",lan->n,pow(lan->tol*lan->quench_tol*lan->quench_tol,(double)(lan->n-1.0)/2.0));
	if(fabs(lan->beta[lan->n])<lan->tol)
	{
		lan->k = -1;
		return -1;
	}
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
	int k;				//the number of converged eigenvalues
	struct zvec* q;			//Lanczos vectors & q[0]=0 q[1] is the initial vector
	struct zvec Ritz;		//Ritz vectors used to orthogonalization
	double* eigs;			//eigenvalues of tridiagonal matrix
	double* eigv;			//eigenvector of tridiagonal matrix
};
*/

void Ritz_i(struct Lanczos* lan, MKL_INT i) // cal the ith Ritz vector
{
	memset(lan->Ritz.vec, 0, (size_t) lan->dim * sizeof(MKL_Complex16));
	for(size_t j=0;j<lan->n;++j)
		z_axpy(&(lan->Ritz), &(lan->q[i+1]), (MKL_Complex16)lan->eigv[i*lan->n+j]);
	return;
}

double quench_tol_Q(struct Lanczos* lan, double dt) 
// cal the tolerance of Lanczos quench iteration for time interval dt
{
	if(lan->n<2)
		return pow(10, 9);
	double min_eig, max_eig, gapdt, tol;
	min_eig = pow(10, 9);
	max_eig = -pow(10,9);
	for(size_t i=0;i<lan->n;++i)
	{
		min_eig = (lan->eigs[i+1]<min_eig)?lan->eigs[i+1]:min_eig;
		max_eig = (lan->eigs[i+1]>max_eig)?lan->eigs[i+1]:max_eig;
	}
	gapdt = 0.25 * dt * (max_eig-min_eig);
	tol = fabs(12.0 * exp(lan->n - gapdt*gapdt/lan->n) * pow(gapdt/lan->n, lan->n));
	tol = (lan->n>2*gapdt)?tol:pow(10,9);
	printf("The %ld-th iteration:\tmin_eig = %e\tmax_eig = %e\t", lan->n, min_eig, max_eig);
	return tol;
}

double first_tol_Q(struct Lanczos* lan, size_t *ind) 
// cal the tolerance of Lanczos quench iteration for time interval dt
{
	double min_eig = lan->eigs[1];
	size_t min_ind = 1;
	for(size_t i=0;i<lan->n;++i)
		min_ind = (lan->eigs[i+1]<min_eig)?(i+1):min_ind;
	min_eig = lan->eigs[min_ind];
	*ind = min_ind;
	double tol = fabs((double)lan->beta[min_ind]*lan->eigv[lan->n*lan->n-lan->n+min_ind-1]);
	printf("The %ld-th iteration:\t",lan->n);
	printf("min_eig = %e\t tol = %e\n",min_eig,tol);
	return tol;
}