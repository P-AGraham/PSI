

bool Lplus(Data_basis *bra_i, size_t j, size_t N_o)// L_j^+
{
    if(j%N_o==N_o-1) return false;
    if(bra_i[j]!=1 || bra_i[j+1]!=0) return false;
    bra_i[j] = 0;
    bra_i[j+1] = 1;
    return true;
}

bool Lminus(Data_basis *bra_i, size_t j, size_t N_o)// L_j^-
{
    if(j%N_o==0) return false;
    if(bra_i[j]!=1 || bra_i[j-1]!=0) return false;
    bra_i[j] = 0;
    bra_i[j-1] = 1;
    return true;
}

bool LplusQ(const Data_basis *bra_i, size_t j, size_t N_o)// L_j^+
{
    if(j%N_o==N_o-1) return false;
    return (bra_i[j]==1) && (bra_i[j+1]==0);
}

bool LminusQ(const Data_basis *bra_i, size_t j, size_t N_o)// L_j^-
{
    if(j%N_o==0) return false;
    return (bra_i[j]==1) && (bra_i[j-1]==0);
}

size_t Lplus(const Data_basis *bra_i, Data_basis *bra_j, size_t j, size_t N_o, size_t N_e)// L_j^+
{
    size_t k = 0;
    for(size_t i=0;i<2*N_o;++i)
    {
        if(bra_i[i])
        {
            bra_j[k] = i;
            if(i==j) bra_j[k] = j+1;
            ++k;
        }
    }
    return dict_num(bra_j, N_e);
}

size_t Lminus(const Data_basis *bra_i, Data_basis *bra_j, size_t j, size_t N_o, size_t N_e)// L_j^-
{
    size_t k = 0;
    for(size_t i=0;i<2*N_o;++i)
    {
        if(bra_i[i])
        {
            bra_j[k] = i;
            if(i==j) bra_j[k] = j-1;
            ++k;
        }
    }
    return dict_num(bra_j, N_e);
}

extern "C" size_t L2Count(MKL_INT *indptr, int num_thread)
{
    indptr[0] = 0;

    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        indptr[i+1] = L2Count_thread(A.basis, A.dim, A.N_o, A.N_e, i);

    for(size_t i=1;i<=A.dim;++i)
        indptr[i] += indptr[i-1];
    return indptr[A.dim];
}

size_t L2Count_thread(Data_basis *k1_basis, size_t k1_dim, \
                    size_t N_o, size_t N_e, size_t i)
{
    size_t count = 0;
    Data_basis *basis_t = k1_basis + i * N_e;
    size_t j1, j2, j3, j4, J12, J34;
    size_t N_orbit = 2 * N_o;
    double s0 = (N_o-1.0)/2;

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    size_t k = 0;

    basis_t = k1_basis + i * N_e;

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
                    ++count;
            }
        }
        /*memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
        for(size_t i=0;i<N_e;++i)
            bra_i[basis_t[i]] = 1;
        if(Lminus(bra_i, u, N_o))
        {
            for(size_t v=0;v<N_orbit;++v)
            {
                if( v!=u && u!=v+1 && LplusQ(bra_i, v, N_o) )
                    ++count;
            }
        }*/
    }

    free(bra_i);
    return count+1;
}

extern "C" size_t L2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread)
{
    indptr[0] = 0;
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
    {
        auto count = L2Generate_thread(indptr, indices, data, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function L2Generate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

size_t L2Generate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, Data_basis *k1_basis, \
                            const hash &hashT, size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    size_t count = 0;
    Data_basis *basis_t = k1_basis + row * N_e;
    size_t j1, j2, j3, j4, J12, J34;
    size_t N_orbit = 2 * N_o;
    double s0 = (N_o-1.0)/2;

    Data_dtype dig = 0.0;
    MKL_INT *indices_thread = indices + indptr[row];
    Data_dtype* data_thread = data + indptr[row];

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    size_t k = 0;

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
                    auto it = hashT.find( bra_num );
                    if(it!=hashT.end())
                    {
                        indices_thread[count] = it->second;
                        //data_thread[count] = 0.5*std::sqrt( (s0+V-s0)*(s0-(V-s0)+1) * (s0-U+s0)*(s0+U-s0+1) ) ;
                        data_thread[count] = std::sqrt( (s0+V-s0)*(s0-(V-s0)+1) * (s0-U+s0)*(s0+U-s0+1) ) ;
                        ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
                        //printf("J-(%lu)J+(%lu): %lu->%lu\n", v, u, row, it->second);
                    }
                }
            }
        }
        /*memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
        for(size_t i=0;i<N_e;++i)
            bra_i[basis_t[i]] = 1;
        if(Lminus(bra_i, u, N_o))
        {
            for(size_t v=0;v<N_orbit;++v)
            {
                if( v!=u && u!=v+1 && LplusQ(bra_i, v, N_o) )
                {
                    size_t V = v%N_o;
                    auto bra_num = Lplus(bra_i, bra_j, v, N_o, N_e);
                    auto it = hashT.find( bra_num );
                    if(it!=hashT.end())
                    {
                        indices_thread[count] = it->second;
                        data_thread[count] = 0.5*std::sqrt( (s0-V+s0)*(s0+(V-s0)+1) * (s0+U-s0)*(s0-(U-s0)+1) ) ;
                        ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
                        //printf("J+(%lu)J-(%lu): %lu->%lu\n", v, u, row, it->second);
                    }
                }
            }
        }*/
    }

    dig = s0*(s0+1)*N_e;
    for(size_t i=0;i<N_e;++i)
    {
        for(size_t j=i+1;j<N_e;++j)
        {
            dig += 2*( (basis_t[i]%N_o)-s0 )*( (basis_t[j]%N_o)-s0 );
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