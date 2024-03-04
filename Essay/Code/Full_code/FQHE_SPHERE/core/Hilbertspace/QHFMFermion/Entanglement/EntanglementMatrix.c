// EntanglementMatrix.c

// | 0,1,...,o_left-1 | o_left,...,N_phi-1 |

int orbitalEntanglementMatrixDimension(MKL_INT *rdmdim, size_t k1, size_t o_left, int num_thread) 
{
	int st;
	if(A.k1 != k1)
	{
		printf("k1(%lu) != A.k1(%lu)", k1, A.k1);
		return -1;
	}

	if(A.partition!=NULL)
		free(A.partition);
	A.partition = (struct bipartition *)malloc(A.dim * sizeof(struct bipartition));
	memset(A.partition, 0, (size_t)A.dim * sizeof(struct bipartition));
	orbitBipartition(A.basis, A.partition, A.dim, A.N_e, o_left, num_thread);
	orbitBipartOrdering(A.basis, A.partition, rdmdim, o_left, num_thread);
	return 0;
}

int orbitalEntanglementMatrix(MKL_Complex16 *ket_k1, MKL_Complex16 *oem, \
								size_t occupy, size_t ky, size_t dim2, int num_thread) 
{
	for (size_t i=0;i<A.dim;++i)
	{
		if(A.partition[i].occupy!=occupy || A.partition[i].ky!=ky)
			continue;
		oem[(size_t)A.partition[i].left*dim2 + A.partition[i].right] += ket_k1[i];
	}
	return 0;
}

int orbitalEntanglementSpetrum(MKL_INT *rdmdim, double *s, MKL_Complex16 *ket_k1, \
								size_t s_length, size_t na, int num_thread) //na=o_left
{
	size_t N_phi = A.N_phi;
	size_t *ind = (size_t *)mkl_malloc((na+2)*N_phi*sizeof(size_t), 64);

	size_t ss = 0;
	for(size_t i=0;i<na+1;++i)
		for(size_t j=0;j<N_phi;++j)
		{
			ind[i*N_phi+j] = ss;
			ss += rdmdim[2*i*N_phi+2*j] * rdmdim[2*i*N_phi+2*j+1];
		}
	MKL_Complex16 *oem = (MKL_Complex16 *)mkl_calloc(ss, sizeof(MKL_Complex16), 64);

	# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<A.dim;++i)
	{
		ky_p = A.partition[i].ky;
		occupy_p = A.partition[i].occupy;
		shift_p = rdmdim[2*occupy_p*N_phi+2*ky_p+1];
		shift_p = A.partition[i].left*shift_p + A.partition[i].right;
		shift_p += ind[occupy_p*N_phi+ky_p];
		oem[shift_p] = ket_k1[i];
	}

	char jobu = 'N';
	char jobv = 'N';
	MKL_Complex16 *u = NULL;
	MKL_Complex16 *vt= NULL;
	size_t m, n;
	int layout = LAPACK_ROW_MAJOR;
	double *superb = (double *)mkl_calloc(ss, sizeof(double), 64);
	double *s_t = s;
	MKL_Complex16 *oem_t = oem;

	for(size_t i=0;i<na+1;++i)
		for(size_t j=0;j<N_phi;++j)
		{
			m = rdmdim[2*i*N_phi+2*j];
			n = rdmdim[2*i*N_phi+2*j+1];
			oem_t = oem + ind[i*N_phi+j];
			LAPACKE_zgesvd(layout, jobu, jobv, m, n, oem_t, n, s_t, u, 1, vt, 1, superb);
			s_t += (m>n)?n:m;
		}

	mkl_free(oem);
	mkl_free(ind);
	mkl_free(superb);
	return 0;
}