


int density(Data_dtype *ket, Data_basis *basis, double *out, \
			size_t N_e, size_t dim)
{
	Data_basis *basis_t = basis;
	double *ket_norm = (double *)malloc((size_t)dim * sizeof(double));

	//# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<dim;++i)
	{
		#ifdef __DATA_DOUBLE
		ket_norm[i] = ket[i]*ket[i];
		#else
		ket_norm[i] = norm(ket[i]);
		#endif
	}

	double temp;
	for(size_t i=0;i<dim;++i)
	{
		temp = ket_norm[i];
		for(size_t j=0;j<N_e;++j)
		{
			out[basis_t[j]] += temp;
		}
		basis_t += N_e;
	}

	free(ket_norm);
	return 0;
}

extern "C" int meausre_density(Data_dtype *ket, double *out)
{
	return density(ket, A.basis, out, A.N_e, A.dim);
}

int sigma_x(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, size_t N_e, size_t N_o, size_t dim)
{
	Data_basis *basis_t = basis;
	Data_dtype *ket_1 = (Data_dtype *)malloc((size_t)dim * sizeof(Data_dtype));

	Data_basis *bra_i = (Data_basis *)malloc(2*N_o * sizeof(Data_basis));
	Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));

	Data_basis *temp = (Data_basis *)calloc(N_o*2, sizeof(unsigned));

	int phase = 0;
	size_t j1, j2;

	for(size_t j=0;j<2*N_o;++j)
	{
		j2 = j;
		j1 = tight_up(j2, N_o);
		memset(ket_1, 0, dim * sizeof(Data_dtype));
		for(size_t i=0;i<dim;++i)
		{
			basis_t = basis + N_e*i;
			memset(temp, 0, N_o*2 * sizeof(Data_basis));
			for(size_t k=0;k<N_e;++k)
				temp[basis_t[k]] = 1;

			//printf("%u->%u: %lu->, temp[j1=%lu]=%u, temp[j2=%lu]=%u\t", j2, j1, i,j1,temp[j1],j2,temp[j2]);
			if(temp[j1]) 
			{
				//printf("kill\n");
				continue;
			}
			//printf("%u->%u: %lu->\n", j2, j1, i);
			if(!temp[j2]) continue;

			auto ind = tight_bond_hopping(basis_t, bra_i, bra_j, N_e, N_o, j1, j2, &phase);
			auto it = hashT.find( ind );
			if(it!=hashT.end())
			{
				ket_1[it->second] += ket[i] * phase;
				//printf("%u->%u: %lu->%lu, phase=%d\n", j2, j1, i, it->second, phase);
				//for(auto ii=0;ii<N_e;++ii)
				//	printf("%u ", bra_j[ii]);
				//printf("\n");
			}
			else
			{
				printf("%lu->%lu: %lu->?, phase=%d\n", j2, j1, i, phase);
				printf("sigma_x: Not found!\n");
				return -1;
			}
		}

		Data_dtype tempc = 0.0;
		for(size_t i=0;i<dim;++i)
		{	
			#ifdef __DATA_DOUBLE
			tempc +=  ket[i]*ket_1[i];
			#else
			tempc +=  conj(ket[i])*ket_1[i];
			#endif
		}
		out[j%N_o] += tempc;
	}

	free(ket_1);
	free(bra_i);
	free(bra_j);
	free(temp);
	return 0;
}

extern "C" int meausre_sigma_x(Data_dtype *ket, Data_dtype *out)
{
	return sigma_x(ket, A.basis, out, A.hashT, A.N_e, A.N_o, A.dim);
}

int Z2(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, size_t N_e, size_t N_o, size_t dim)
{
	Data_basis *basis_t = basis;
	Data_dtype *ket_1 = (Data_dtype *)malloc((size_t)dim * sizeof(Data_dtype));
	memset(ket_1, 0, dim * sizeof(Data_dtype));
	char phase = 1;

	for(size_t i=0;i<dim;++i)
	{
		auto it = hashT.find( base_translation(basis_t, N_o, 2*N_o, N_e, &phase) );
		if(it!=hashT.end())
		{
			ket_1[it->second] = ket[i] * phase;
			//printf("%lu->%lu: phase=%d\n", i, it->second, phase);
		}
		else
		{
			return -1;
		}
		basis_t += N_e;
	}

	Data_dtype tempc = 0.0;
	for(size_t i=0;i<dim;++i)
	{
		//std::cout << ket[i] << "\t" << ket_1[i] << std::endl;
		#ifdef __DATA_DOUBLE
		tempc +=  ket[i]*ket_1[i];
		#else
		tempc +=  conj(ket[i])*ket_1[i];
		#endif
	}
	out[0] = tempc;

	free(ket_1);
	return 0;
}


extern "C" int meausre_Z2(Data_dtype *ket, Data_dtype *out)
{
	return Z2(ket, A.basis, out, A.hashT, A.N_e, A.N_o, A.dim);
}

extern "C" int meausre_PH(Data_dtype *ket, Data_dtype *out)
{
	return PH(ket, A.basis, out, A.hashT, A.N_e, A.N_o, A.dim);
}

int PH(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, size_t N_e, size_t N_o, size_t dim)
{
	Data_basis *basis_t = basis;
	Data_dtype *ket_1 = (Data_dtype *)malloc((size_t)dim * sizeof(Data_dtype));
	Data_basis *PHbasis = (Data_basis *)malloc((size_t)N_e * sizeof(Data_basis));
	memset(ket_1, 0, dim * sizeof(Data_dtype));
	char phase = 1;

	for(size_t i=0;i<dim;++i)
	{
		PH_transformation(basis_t, PHbasis, N_o, N_e, &phase);
		//print_basis_mod(basis_t, N_e, N_o);
		//printf("\t->\t");
		//print_basis_mod(PHbasis, N_e, N_o);
		auto it = hashT.find( dict_num(PHbasis, 2*N_o-N_e) );
		if(it!=hashT.end())
		{
			ket_1[it->second] = ket[i] * phase;
			//printf("%lu->%lu: phase=%d\n", i, it->second, phase);
		}
		else
		{
			return -1;
		}
		basis_t += N_e;
	}

	Data_dtype tempc = 0.0;
	for(size_t i=0;i<dim;++i)
	{
		//std::cout << ket[i] << "\t" << ket_1[i] << std::endl;
		#ifdef __DATA_DOUBLE
		tempc +=  ket[i]*ket_1[i];
		#else
		tempc +=  conj(ket[i])*ket_1[i];
		#endif
	}
	out[0] = tempc;

	free(PHbasis);
	free(ket_1);
	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<Z2>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" int meausre_PH_inZ2space(Data_dtype *ket, Data_dtype *out)
{
	return PH(ket, A.basis, out, A.hashT, A.Z2_phase, A.Z_2, A.N_e, A.N_o, A.dim);
}

int PH(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, char* Z2_phase, int Z_2, size_t N_e, size_t N_o, size_t dim)
{
	Data_basis *basis_t = basis;
	Data_dtype *ket_1 = (Data_dtype *)malloc((size_t)dim * sizeof(Data_dtype));
	Data_basis *PHbasis = (Data_basis *)malloc((size_t)N_e * sizeof(Data_basis));
	memset(ket_1, 0, dim * sizeof(Data_dtype));
	char phase = 1;
	size_t ind;
	int n_find;

	for(size_t i=0;i<dim;++i)
	{
		PH_transformation(basis_t, PHbasis, N_o, N_e, &phase);
		//auto it = hashT.find( dict_num(PHbasis, 2*N_o-N_e) );
		ind = find_Z2_index(PHbasis, hashT, N_e, N_o, n_find);
        if(ind!=ULIMAX)
		{
            if( n_find==2 && Z_2*Z2_phase[ind]<0 )
                phase *= -1;
			ket_1[ind] = ket[i] * std::sqrt( std::abs((double)Z2_phase[i]/Z2_phase[ind]) ) * phase;
			//printf("%lu->%lu: phase=%d\n", i, it->second, phase);
		}
		basis_t += N_e;
	}

	Data_dtype tempc = 0.0;
	for(size_t i=0;i<dim;++i)
	{
		//std::cout << ket[i] << "\t" << ket_1[i] << std::endl;
		#ifdef __DATA_DOUBLE
		tempc +=  ket[i]*ket_1[i];
		#else
		tempc +=  conj(ket[i])*ket_1[i];
		#endif
	}
	out[0] = tempc;

	free(PHbasis);
	free(ket_1);
	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<Z2PH>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
