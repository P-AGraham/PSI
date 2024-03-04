
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
	int k;					//the number of converged eigenvalues
	struct zvec* q;			//Lanczos vectors & q[0]=0 q[1] is the initial vector
	MKL_INT qn;				//the number of allocated q vector
	struct zvec Ritz;		//Ritz vectors used to orthogonalization
	double* eigs;			//eigenvalues of tridiagonal matrix
	double* eigv;			//eigenvector of tridiagonal matrix
	double quench_tol;		//the tolerance of quench time
};

struct Lanczos create_lanczos(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
							MKL_INT dim, double tol, MKL_Complex16* initv, MKL_INT Hermitian);
							// 0 General Matrix 1 Hermitian matrix

void destroy_lanczos(struct Lanczos* lan);

int Lanczos_realloc(struct Lanczos* lan);

//---------------------------LanczosReO/LanczosReothogonalization.c----------------------------------------------------

int LanczosIteration(struct Lanczos* lan); 	// update the Lanczos struct lan

void Ritz_i(struct Lanczos* lan, MKL_INT i); // cal the ith Ritz vector

double quench_tol_Q(struct Lanczos* lan, double dt) ;
// cal the tolerance of Lanczos quench iteration for time interval dt

double first_tol_Q(struct Lanczos* lan, size_t *ind);
// cal the lowest eigenvalue/vector's error bound