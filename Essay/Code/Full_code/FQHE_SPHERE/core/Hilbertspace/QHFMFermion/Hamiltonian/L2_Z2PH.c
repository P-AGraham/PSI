
extern "C" size_t L2Z2PHCount(MKL_INT *indptr, int num_thread)
{
    indptr[0] = 0;

    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        indptr[i+1] = L2Z2PHCount_thread(A.basis, A.PHindex, A.hashT, A.PHhashT, A.Z2_phase, A.PH_phase,\
                                         A.Z2PH_length, A.Z_2, A.PH, A.dim, A.N_o, A.N_e, i);

    for(size_t i=1;i<=A.dim;++i)
        indptr[i] += indptr[i-1];
    return indptr[A.dim];
}

size_t L2Z2PHCount_thread(Data_basis *k1_basis, const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                             const char* Z2_phase, const char* PH_phase, const char* Z2PH_length, int Z_2, int PH, \
                             size_t k1_dim, size_t N_o, size_t N_e, size_t i)
{
    size_t count = 0;
    const Data_basis *basis_t = k1_basis + PHindex[i] * N_e;
    size_t j1, j2, j3, j4, J12, J34, ind;
    size_t N_orbit = 2 * N_o;
    char cphase;
    int iphase;
    int length;
    double s0 = (N_o-1.0)/2;

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    size_t k = 0;

    basis_t = k1_basis + PHindex[i] * N_e;

    for(size_t u=0;u<N_orbit;++u)
    {
        memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
        for(size_t i=0;i<N_e;++i)
            bra_i[basis_t[i]] = 1;
        if(Lplus(bra_i, u, N_o))
        {
            for(size_t v=0;v<N_orbit;++v)
            {
                if( v!=u && v!=u+1 &&LminusQ(bra_i, v, N_o) )
                {
                    auto bra_num = Lminus(bra_i, bra_j, v, N_o, N_e);
                    //ind = find_Z2_index(bra_j, hashT, N_e, N_o, n_find);
                    ind = find_Z2PH_index(bra_j, hashT, PHhashT, Z2_phase, PH_phase, Z2PH_length,\
                                        PHindex, N_e, N_o, length, cphase);
                    if(ind!=ULIMAX)
                        ++count;
                }
            }
        }
    }

    free(bra_i);
    free(bra_j);
    return count+1;
}

extern "C" size_t L2Z2PHGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread)
{
    indptr[0] = 0;
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
    {
        auto count = L2Z2PHGenerate_thread(indptr, indices, data, A.basis, A.PHindex, A.hashT, A.PHhashT,\
                                         A.Z2_phase, A.PH_phase, A.Z2PH_length, A.Z_2, A.PH, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function L2Z2Generate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

//#define __HANCHAO

size_t L2Z2PHGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    Data_basis *k1_basis, const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                    const char* Z2_phase, const char* PH_phase, const char* Z2PH_length, int Z_2, int PH, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    size_t count = 0;
    const Data_basis *basis_t = k1_basis + PHindex[row] * N_e;
    size_t j1, j2, j3, j4, J12, J34, ind;
    size_t N_orbit = 2 * N_o;
    int bralength, ketlength;
    char cphase;
    double s0 = (N_o-1.0)/2;

    Data_dtype dig = 0.0;
    MKL_INT *indices_thread = indices + indptr[row];
    Data_dtype* data_thread = data + indptr[row];

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    size_t k = 0;

    #ifdef __HANCHAO
    print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    find_Z2PH_index(basis_t, hashT, PHhashT, Z2_phase, PH_phase, Z2PH_length, PHindex, N_e, N_o, ketlength, cphase);
    for(size_t u=0;u<N_orbit;++u)
    {
        memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
        for(size_t i=0;i<N_e;++i)
            bra_i[basis_t[i]] = 1;
        size_t U = u%N_o;
        if(Lplus(bra_i, u, N_o))
        {
            for(size_t v=0;v<N_orbit;++v)
            {
                if( v!=u && v!=u+1 &&LminusQ(bra_i, v, N_o) )
                {
                    size_t V = v%N_o;
                    auto bra_num = Lminus(bra_i, bra_j, v, N_o, N_e);
                    //ind = find_Z2_index(bra_j, hashT, N_e, N_o, n_find);
                    ind = find_Z2PH_index(bra_j, hashT, PHhashT, Z2_phase, PH_phase, Z2PH_length,\
                                            PHindex, Z_2, PH, N_e, N_o, bralength, cphase);
                    if(ind!=ULIMAX)
                    {
                        indices_thread[count] = ind;
                        //data_thread[count] = 0.5*std::sqrt( (s0+V-s0)*(s0-(V-s0)+1) * (s0-U+s0)*(s0+U-s0+1) ) ;
                        data_thread[count] = std::sqrt( (s0+V-s0)*(s0-(V-s0)+1) * (s0-U+s0)*(s0+U-s0+1) ) ;
                        //if( n_find==2 && Z_2*Z2_phase[ind]<0 )
                        //    data_thread[count] *= -1;
                        data_thread[count] *= std::sqrt( std::abs( (double)ketlength/bralength ) )*cphase;
                        ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
                        //printf("J-(%lu)J+(%lu): %lu->%lu\n", v, u, row, it->second);
                    }
                }
            }
        }
    }

    dig = s0*(s0+1)*N_e;
    double digg = 0;
    for(size_t i=0;i<N_e;++i)
    {
        for(size_t j=i+1;j<N_e;++j)
        {
            dig += 2*( (basis_t[i]%N_o)-s0 )*( (basis_t[j]%N_o)-s0 );
            #ifdef __HANCHAO
            digg += 2*( (basis_t[i]%N_o)-s0 )*( (basis_t[j]%N_o)-s0 );
            printf("\t2*L^z_%lu*L^z_%lu = 2*%lf*%lf = %lf, total = %lf\n", basis_t[i], basis_t[j], (basis_t[i]%N_o)-s0, (basis_t[j]%N_o)-s0, 2*( (basis_t[i]%N_o)-s0 )*( (basis_t[j]%N_o)-s0 ), digg );
            #endif
        }
    }
    for(size_t i=0;i<N_e-1;++i)
    {
        if(basis_t[i]%N_o==N_o-1) continue;
        size_t U = basis_t[i]%N_o;
        if(basis_t[i]+1==basis_t[i+1])
        {
            dig -= (s0-U+s0)*(s0+U-s0+1);
        }
    }
    indices_thread[count] = row;
    data_thread[count] = dig;

    free(bra_i);
    free(bra_j);
    return count+1;
}