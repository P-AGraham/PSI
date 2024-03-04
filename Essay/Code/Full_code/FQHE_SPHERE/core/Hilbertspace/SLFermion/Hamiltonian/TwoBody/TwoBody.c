// FermionHamiltonianSymmetry.c

size_t twoBodyCount(MKL_INT *indptr, int num_thread)
{
	indptr[0] = 0;

	# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<A.dim;++i)
		indptr[i+1] = twoBodyCount_thread(i);

	for(size_t i=1;i<=A.dim;++i)
		indptr[i] += indptr[i-1];
	return indptr[A.dim];
}

size_t twoBodyCount_thread(size_t row)
{
	size_t bra_ind, bra_num, N_phi, N_e, n_data = 0;
	int j1, j2, j3, j4, J12, J34, phase;
	N_phi = A.N_o;
	N_e = A.N_e;

	MKL_INT *index = A.index;
	unsigned char *basis_t = A.basis + row * N_e;

	unsigned char *bra_i = (unsigned char *)malloc(N_phi * sizeof(unsigned char));
	unsigned char *bra_j = (unsigned char *)malloc(N_e * sizeof(unsigned char));
	unsigned char *diff = (unsigned char *)malloc((N_phi-N_e) * sizeof(unsigned char));
	size_t k = 0;
	for(size_t j=0;j<N_phi;++j)
	{
		if(*basis_t == j)
			++basis_t;
		else
		{
			diff[k] = j;
			++k;
		}
	}
	basis_t = A.basis + row * N_e;

	for(size_t u=0;u<N_e-1;++u)
	{
		j3 = basis_t[u];
		for(size_t v=u+1;v<N_e;++v)
		{
			j4 = basis_t[v];
			J34 = j3 + j4;
			for(size_t s=0;s<N_phi-N_e-1;++s)
			{
				j2 = diff[s];
				for(size_t t=s+1;t<N_phi-N_e;++t)
				{
					j1 = diff[t];
					J12 = j1 + j2;
					if(J12!=J34)
						continue;
					bra_num = hopping2Body(basis_t, bra_i, bra_j, N_phi, N_e, j1, j2, j3, j4, &phase);
					bra_ind = index[bra_num];//find_node(head, bra_num, hash_size);
					if((A.HermitianQ!=0)&&(bra_ind<row))
						continue;
					//printf("%u\n",bra_ind);
					++n_data;
				}
			}
		}
	}
	++n_data;

	free(bra_i);
	free(bra_j);
	free(diff);
	return n_data;
}

size_t twoBodyGenerate(double *A_list, MKL_INT *indptr, MKL_INT *indices, \
						Data_dtype* data, int num_thread)
{
	indptr[0] = 0;
	//# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<A.dim;++i)
	{
		if(twoBodyGenerate_thread(i, A_list, indptr, indices, data)!= indptr[i+1] - indptr[i])
			exit(10);
	}

	return indptr[A.dim];
}

size_t twoBodyGenerate_thread(size_t row, double *A_list, MKL_INT *indptr, \
								MKL_INT *indices, Data_dtype *data)
{
	size_t bra_ind, bra_num, N_phi, N_e, u, n_data = 0;
	int j1, j2, j3, j4, J12, J34, phase;
	N_phi = A.N_o;
	N_e = A.N_e;
	u = N_phi;

	Data_dtype dig = 0.0;
	MKL_INT *index = A.index;
	unsigned char *basis_t = A.basis + row * N_e;
	MKL_INT *indices_thread = indices + indptr[row];
	Data_dtype *data_thread = data + indptr[row];

	unsigned char *bra_i = (unsigned char *)malloc(N_phi * sizeof(unsigned char));
	unsigned char *bra_j = (unsigned char *)malloc(N_e * sizeof(unsigned char));
	unsigned char *diff = (unsigned char *)malloc((N_phi-N_e) * sizeof(unsigned char));
	size_t k = 0;
	for(size_t j=0;j<N_phi;++j)
	{
		if(*basis_t == j)
			++basis_t;
		else
		{
			diff[k] = j;
			++k;
		}
	}
	basis_t = A.basis + row * N_e;

	//for(size_t i=0;i<N_e;++i)
	//	printf("%d\t", basis_t[i]);
	//printf("\n");

	for(size_t u1=0;u1<N_e-1;++u1)
	{
		j3 = basis_t[u1];
		for(size_t v=u1+1;v<N_e;++v)
		{
			j4 = basis_t[v];
			J34 = j3 + j4;
			for(size_t s=0;s<N_phi-N_e-1;++s)
			{
				j2 = diff[s];
				for(size_t t=s+1;t<N_phi-N_e;++t)
				{
					j1 = diff[t];
					J12 = j1 + j2;
					if(J12!=J34)
						continue;
					bra_num = hopping2Body(basis_t, bra_i, bra_j, N_phi, N_e, j1, j2, j3, j4, &phase);
					bra_ind = index[bra_num];//find_node(head, bra_num, hash_size);
					if((A.HermitianQ!=0)&&(bra_ind<row))
						continue;
					dig = 2.0 * A_list[j1*u*u*u + j2*u*u + j3*u + j4];
					dig -=2.0 * A_list[j2*u*u*u + j1*u*u + j3*u + j4];
					indices_thread[n_data] = bra_ind;
					data_thread[n_data] = dig * phase;
					//printf("\t%d %d %d %d:\t%lu %lf %lf %lf %lf %lf %lf %d\n", j1, j2, j3, j4, bra_ind,\
						A_list[j1*u*u*u + j2*u*u + j3*u + j4],\
						A_list[j2*u*u*u + j1*u*u + j4*u + j3],\
						A_list[j1*u*u*u + j2*u*u + j4*u + j3],\
						A_list[j2*u*u*u + j1*u*u + j3*u + j4],\
						A_list[7034],\
						dig,phase);
					++n_data;
				}
			}
		}
	}
	dig = 0.0;
	for(size_t s1=0;s1<N_e-1;++s1)
	{
		for(size_t s2=s1+1;s2<N_e;++s2)
		{
			j3 = basis_t[s2];
			j4 = basis_t[s1];
			j1 = j4;
			j2 = j3;
			dig += 2.0 * A_list[j1*u*u*u + j2*u*u + j3*u + j4];
			dig -= 2.0 * A_list[j2*u*u*u + j1*u*u + j3*u + j4];
		}
	}
	indices_thread[n_data] = row;
	data_thread[n_data] = dig;
	++n_data;

	free(bra_i);
	free(bra_j);
	free(diff);
	return n_data;
}

size_t hopping2Body(unsigned char *basis, unsigned char *bra_i, unsigned char *bra_j, \
					size_t N_phi, size_t N_e, int j1, int j2, int j3, int j4, int *phase)
{
	// j4 > j3, j1 > j2

	// initialization
	memset(bra_i, 0, (size_t)N_phi * sizeof(unsigned char));

	// translation occupy -> binary
	for(size_t i=0;i<N_e;++i)
		++bra_i[basis[i]];
	*phase = 0;

	// annihilation
	for(size_t i=j3+1;i<j4;++i)
		*phase += bra_i[i];
	--bra_i[j4];
	--bra_i[j3];

	// creation
	for(size_t i=j2+1;i<j1;++i)
		*phase += bra_i[i];
	++bra_i[j2];
	++bra_i[j1];

	*phase = (*phase&1)?-1:1;
	
	// translation binary -> occupy
	size_t k = 0;
	for(size_t i=0;i<N_phi;++i)
	{
		if(bra_i[i])
		{
			bra_j[k] = i;
			++k;
		}
	}

	return dict_num(bra_j, N_e);
}