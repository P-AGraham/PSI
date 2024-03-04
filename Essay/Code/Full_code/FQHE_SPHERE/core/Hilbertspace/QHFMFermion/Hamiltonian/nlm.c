
extern "C" size_t n00ACount(MKL_INT *indptr, Data_dtype *A_matirx, int num_thread)
{
    indptr[0] = 0;

    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        indptr[i+1] = n00ACount_thread(A_matirx, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);

    for(size_t i=1;i<=A.dim;++i)
        indptr[i] += indptr[i-1];
    return indptr[A.dim];
}

//#define __DEBUG_N00A_COUNT

size_t n00ACount_thread(const Data_dtype *A_matirx, Data_basis *basis, const hash & hashT, size_t k1_dim, size_t N_o, size_t N_e, size_t i)
{
    size_t count = 0;
    const Data_basis *basis_t = basis + i * N_e;
    size_t ind, N_orbit = 2 * N_o;
    char phase;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
    
    #ifdef __DEBUG_N00A_COUNT
    print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    if( std::abs(A_matirx[1])>1e-15 )
    {
        for(size_t m=0;m<N_o;++m)
        {
            if(hopping_DowntoUp(bin_basis, bra_basis, m, N_e, N_o, phase))
            {
                ind = dict_num(bra_basis, N_e);
                auto it = hashT.find(ind);
                if(it!=hashT.end())
                    ++count;
                #ifdef __DEBUG_N00A_COUNT
                printf("\t");
                printf("hopping_DowntoUp %lu:\t", m);
                print_basis(bra_basis, N_e, N_o);
                printf("\t resu ind: %lu\n", ind);
                #endif
            }
        }
    }

    if( std::abs(A_matirx[2])>1e-15 )
    {
        for(size_t m=0;m<N_o;++m)
        {
            if(hopping_UptoDown(bin_basis, bra_basis, m, N_e, N_o, phase))
            {
                ind = dict_num(bra_basis, N_e);
                auto it = hashT.find(ind);
                if(it!=hashT.end())
                    ++count;
                #ifdef __DEBUG_N00A_COUNT
                printf("\t");
                printf("hopping_UptoDown %lu:\t", m);
                print_basis(bra_basis, N_e, N_o);
                printf("\t resu ind: %lu\n", ind);
                #endif
            }
        }
    }

    if( std::abs(A_matirx[0])>1e-15 || std::abs(A_matirx[3])>1e-15 )
        count += 1;

    free(bin_basis);
    free(bra_basis);
    return count;
}

extern "C" size_t n00AGenerate(Data_dtype *A_matirx, MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread)
{
    indptr[0] = 0;
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
    {
        auto count = n00AGenerate_thread(A_matirx, indptr, indices, data, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function n00AGenerate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

//#define __DEBUG_N00A_GENERATE

size_t n00AGenerate_thread(const Data_dtype *A_matirx, MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                            Data_basis *basis, const hash & hashT, size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    size_t count = 0;
    const Data_basis *basis_t = basis + row * N_e;
    size_t ind, N_orbit = 2 * N_o;
    char phase, cphase;
    int length;
    double comm = 1/std::sqrt(4*M_PI);
    Data_dtype dig = 0;

    MKL_INT *indices_thread = indices + indptr[row];
    Data_dtype* data_thread = data + indptr[row];

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
    
    #ifdef __DEBUG_N00A_GENERATE
    print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    if( std::abs(A_matirx[1])>1e-15 )
    {
        for(size_t m=0;m<N_o;++m)
        {
            if(hopping_DowntoUp(bin_basis, bra_basis, m, N_e, N_o, cphase))
            {
                ind = dict_num(bra_basis, N_e);
                auto it = hashT.find(ind);
                if(it!=hashT.end())
                {
                    phase = cphase;
                    ind = it->second;
                    //phase *= (m&1)?-1:1;  // int code m = m + s0 in equation
                    indices_thread[count] = ind;
                    data_thread[count] = A_matirx[1] * phase * comm;
                    ++count;
                    #ifdef __DEBUG_N00A_GENERATE
                    printf("\t");
                    printf("hopping_DowntoUp %lu->%lu:\t", m+N_o, m);
                    print_basis(bra_basis, N_e, N_o);
                    printf("\t resu ind: %lu, phase=%+d, resu=%lf\n", ind, phase, data_thread[count-1]);
                    #endif
                }
            }
        }
    }

    if( std::abs(A_matirx[2])>1e-15 )
    {
        for(size_t m=0;m<N_o;++m)
        {
            if(hopping_UptoDown(bin_basis, bra_basis, m, N_e, N_o, cphase))
            {
                ind = dict_num(bra_basis, N_e);
                auto it = hashT.find(ind);
                if(it!=hashT.end())
                {
                    phase = cphase;
                    ind = it->second;
                    //phase *= (m&1)?-1:1;  // int code m = m + s0 in equation
                    indices_thread[count] = ind;
                    data_thread[count] = A_matirx[2] * phase * comm;
                    ++count;
                    #ifdef __DEBUG_N00A_GENERATE
                    printf("\t");
                    printf("hopping_UptoDown %lu->%lu:\t", m, m+N_o);
                    print_basis(bra_basis, N_e, N_o);
                    printf("\t resu ind: %lu, phase=%+d, resu=%lf\n", ind, phase, data_thread[count-1]);
                    #endif
                }
            }
        }
    }

    dig = 0.0;
    if( std::abs(A_matirx[0])>1e-15 || std::abs(A_matirx[3])>1e-15 )
    {
        for(size_t m=0;m<N_o;++m)
        {
            dig += bin_basis[m] * A_matirx[0];
            dig += bin_basis[m+N_o] * A_matirx[3];
        }
        indices_thread[count] = row;
        data_thread[count] = dig * comm;
        count += 1;
    }

    free(bin_basis);
    free(bra_basis);
    return count;
}


// basis:
//          up:     <N_o
//          down:   >=N_o

// c^dagger_down c_up
bool hopping_UptoDown(Data_basis* bin_basis, Data_basis* resu_basis, size_t m, size_t N_e, size_t N_o, char& phase)
{
    size_t N_orbit = 2*N_o;
    bool resu = false;
    phase = 0;
    if(bin_basis[m+N_o]==0&&bin_basis[m]==1)
    {
        bin_basis[m] = 0;
        bin_basis[m+N_o] = 1;

        resu = true;
        for(auto i=m+1;i<m+N_o;++i)
            phase += bin_basis[i];
        phase = (phase&1)?-1:1;
        binary_to_basis(bin_basis, resu_basis, N_orbit);

        bin_basis[m] = 1;
        bin_basis[m+N_o] = 0;
    }

    return resu;
}

bool hopping_DowntoUp(Data_basis* bin_basis, Data_basis* resu_basis, size_t m, size_t N_e, size_t N_o, char& phase)
{
    size_t N_orbit = 2*N_o;
    bool resu = false;
    phase = 0;
    if(bin_basis[m]==0&&bin_basis[m+N_o]==1)
    {
        bin_basis[m] = 1;
        bin_basis[m+N_o] = 0;

        resu = true;
        for(auto i=m+1;i<m+N_o;++i)
            phase += bin_basis[i];
        phase = (phase&1)?-1:1;
        binary_to_basis(bin_basis, resu_basis, N_orbit);

        bin_basis[m] = 0;
        bin_basis[m+N_o] = 1;
    }

    return resu;
}

// c^dagger_m1 c_{m1-m}
bool hopping_UptoDown(Data_basis* bin_basis, Data_basis* resu_basis, int m1, int m, size_t N_e, size_t N_o, char& phase)
{
    phase = 1;
    if( m1-m<0 || m1-m>=(int)N_o )
        return false;
    size_t N_orbit = 2*N_o;
    bool resu = false;
    phase = 0;
    if(bin_basis[m1+N_o]==0&&bin_basis[m1-m]==1)
    {
        bin_basis[m1-m] = 0;
        bin_basis[m1+N_o] = 1;

        resu = true;
        for(auto i=m1-m+1;i<m1+N_o;++i)
            phase += bin_basis[i];
        phase = (phase&1)?-1:1;
        binary_to_basis(bin_basis, resu_basis, N_orbit);

        bin_basis[m1-m] = 1;
        bin_basis[m1+N_o] = 0;
    }

    return resu;
}

bool hopping_DowntoUp(Data_basis* bin_basis, Data_basis* resu_basis, int m1, int m, size_t N_e, size_t N_o, char& phase)
{
    phase = 1;
    if( m1-m<0 || m1-m>=(int)N_o )
        return false;
    phase = 0;
    size_t N_orbit = 2*N_o;
    bool resu = false;
    if(bin_basis[m1+N_o]==1&&bin_basis[m1-m]==0)
    {
        bin_basis[m1-m] = 1;
        bin_basis[m1+N_o] = 0;

        resu = true;
        for(auto i=m1-m+1;i<m1+N_o;++i)
            phase += bin_basis[i];
        phase = (phase&1)?-1:1;
        binary_to_basis(bin_basis, resu_basis, N_orbit);

        bin_basis[m1-m] = 0;
        bin_basis[m1+N_o] = 1;
    }

    return resu;
}

/*
    layerm1      layerm    hopping_type            hopping                        hoppping_in_code
        0           0       Up_to_Up       c^dagger_{up, m1}c_{up,m1-m}        c^dagger_{m1}c_{m1-m}
        1           0       Up_to_Dn       c^dagger_{dn, m1}c_{up,m1-m}      c^dagger_{m1+N_o}c_{m1-m}
        0           1       Dn_to_Up       c^dagger_{up, m1}c_{dn,m1-m}      c^dagger_{m1}c_{m1-m+N_o}
        1           1       Dn_to_Dn       c^dagger_{dn, m1}c_{dn,m1-m}     c^dagger_{m1+N_o}c_{m1-m+N_o}
*/
bool hopping_LayertoLayer(Data_basis* bin_basis, Data_basis* resu_basis, int m1, int m, int layerm1, int layerm,
                             size_t N_e, size_t N_o, char& phase)
{
    phase = 1;
    if( m1-m<0 || m1-m>=(int)N_o )
        return false;
    size_t N_orbit = 2*N_o;
    bool resu = false;
    phase = 0;
    if(bin_basis[m1+N_o*layerm1]==1&&bin_basis[m1-m+N_o*layerm]==0)
    {
        bin_basis[m1-m+N_o*layerm] = 1;
        bin_basis[m1+N_o*layerm1] = 0;

        resu = true;
        for(auto i=std::min(m1-m+N_o*layerm, m1+N_o*layerm1)+1;i<std::max(m1-m+N_o*layerm, m1+N_o*layerm1);++i)
            phase += bin_basis[i];
        phase = (phase&1)?-1:1;
        binary_to_basis(bin_basis, resu_basis, N_orbit);

        bin_basis[m1-m+N_o*layerm] = 0;
        bin_basis[m1+N_o*layerm1] = 1;
    }

    return resu;
}
bool hopping_LayertoLayer(Data_basis* bin_basis, int m1, int m, int layerm1, int layerm, size_t N_e, size_t N_o, char& phase)
{
    phase = 1;
    if( m1-m<0 || m1-m>=(int)N_o )
        return false;
    if( m==0 && layerm1==layerm )
    {
        phase = 1;
        return bin_basis[m1+N_o*layerm1]==1;
    }
    phase = 0;

    size_t N_orbit = 2*N_o;
    bool resu = false;
    if(bin_basis[m1+N_o*layerm1]==1&&bin_basis[m1-m+N_o*layerm]==0)
    {
        bin_basis[m1-m+N_o*layerm] = 1;
        bin_basis[m1+N_o*layerm1] = 0;

        resu = true;
        for(auto i=std::min(m1-m+N_o*layerm, m1+N_o*layerm1)+1;i<std::max(m1-m+N_o*layerm, m1+N_o*layerm1);++i)
            phase += bin_basis[i];
        phase = (phase&1)?-1:1;
    }

    return resu;
}

void restore_LayertoLayer(Data_basis* bin_basis, int m1, int m, int layerm1, int layerm, size_t N_o)
{
    if( m1-m<0 || m1-m>=(int)N_o )
        return;
    if( m==0 && layerm1==layerm )
        return;

    if(bin_basis[m1+N_o*layerm1]==0&&bin_basis[m1-m+N_o*layerm]==1)
    {
        bin_basis[m1-m+N_o*layerm] = 0;
        bin_basis[m1+N_o*layerm1] = 1;
    }
    return;
}

extern "C" size_t D00ACount(MKL_INT *indptr, Data_dtype *A_matirx, Data_dtype *B_matirx, int num_thread)
{
    indptr[0] = 0;

    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        indptr[i+1] = D00ACount_thread(A_matirx, B_matirx, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);

    for(size_t i=1;i<=A.dim;++i)
        indptr[i] += indptr[i-1];
    return indptr[A.dim];
}

//#define __DEBUG_D00A_COUNT

size_t D00ACount_thread(const Data_dtype *A_matirx, const Data_dtype *B_matirx, Data_basis *basis, const hash & hashT,
                         size_t k1_dim, size_t N_o, size_t N_e, size_t i)
{
    int layers = 2;
    size_t count = 0;
    const Data_basis *basis_t = basis + i * N_e;
    size_t ind, N_orbit = layers * N_o;
    char phase;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
    
    #ifdef __DEBUG_D00A_COUNT
    print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    for(int m=-(int)N_o;m<=(int)N_o;++m)
    {
        for(int j1=0;j1<layers*layers;j1++)
        {
            for(int j2=0;j2<layers*layers;j2++)
            {
                if( std::abs(A_matirx[j1])>1e-15 && std::abs(B_matirx[j2])>1e-15 )
                {
                    int layerm1 = j1/2;
                    int layerm1m = j1&1;
                    int layerm2 = j2/2;
                    int layerm2m = j2&1;
                    for(int m1=0;m1<(int)N_o;++m1)
                    {
                        for(int m2=0;m2<(int)N_o;++m2)
                        {
                            basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
                            if(!hopping_LayertoLayer(bin_basis, m1, -m, layerm1, layerm1m, N_e, N_o, phase))
                                continue;
                            if(!hopping_LayertoLayer(bin_basis, m2, m, layerm2, layerm2m, N_e, N_o, phase))
                                continue;

                            binary_to_basis(bin_basis, bra_basis, N_orbit);

                            ind = dict_num(bra_basis, N_e);
                            auto it = hashT.find(ind);
                            if(it!=hashT.end())
                            {
                                ++count;
                            }
                            #ifdef __DEBUG_D00A_COUNT
                            //if(i!=0)
                            //    continue;
                            printf("\thopping");
                            printf(" %d->%d ", m1+layerm1*N_o, m1+m+layerm1m*N_o);
                            printf("and %d->%d(m1=%d, m2=%d, m=%d):\t", m2+layerm2*N_o, m2-m+layerm2m*N_o, m1, m2,m);
                            print_basis(bra_basis, N_e, N_o);
                            printf("(j1=%d,j2=%d)(layerm1=%d,layerm1m=%d,layerm2=%d,layerm2m=%d)",j1,j2,layerm1,layerm1m,layerm2,layerm2m);
                            if(it!=hashT.end())
                                printf("\t resu num: %lu, ind: %lu\n", ind, it->second);
                            else
                                printf("\t resu num: %lu, ind: hashT.end()\n", ind);
                            #endif
                        }
                    }
                }
            }
        }
    }

    free(bin_basis);
    free(bra_basis);
    return count;
}

extern "C" size_t D00AGenerate(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread)
{
    indptr[0] = 0;
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
    {
        auto count = D00AGenerate_thread(A_matirx, B_matirx, A_list, indptr, indices, data, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function D00AGenerate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

//#define __DEBUG_D00A_GENERATE

size_t D00AGenerate_thread(const Data_dtype *A_matirx, const Data_dtype *B_matirx, const Data_dtype *A_list, 
                            MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, Data_basis *basis, const hash & hashT, 
                            size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    int layers = 2;
    size_t count = 0;
    const Data_basis *basis_t = basis + row * N_e;
    size_t ind, N_orbit = layers * N_o;
    char phase, cphase;
    int length;
    //double comm = 1/std::sqrt(4*M_PI);
    //Data_dtype dig = 0;

    MKL_INT *indices_thread = indices + indptr[row];
    Data_dtype* data_thread = data + indptr[row];

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
    
    #ifdef __DEBUG_D00A_GENERATE
    printf("row=%lu\t", row);print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    for(int m=1-(int)N_o;m<(int)N_o;++m)
    {
        for(int j1=0;j1<layers*layers;j1++)
        {
            for(int j2=0;j2<layers*layers;j2++)
            {
                if( std::abs(A_matirx[j1])>1e-15 && std::abs(B_matirx[j2])>1e-15 )
                {
                    int layerm1 = j1/2;
                    int layerm1m = j1&1;
                    int layerm2 = j2/2;
                    int layerm2m = j2&1;
                    for(int m1=0;m1<(int)N_o;++m1)
                    {
                        for(int m2=0;m2<(int)N_o;++m2)
                        {
                            basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
                            if(!hopping_LayertoLayer(bin_basis, m1, -m, layerm1, layerm1m, N_e, N_o, cphase))
                                continue;
                            phase = cphase;
                            if(!hopping_LayertoLayer(bin_basis, m2, m, layerm2, layerm2m, N_e, N_o, phase))
                                continue;
                            phase *= cphase;

                            binary_to_basis(bin_basis, bra_basis, N_orbit);

                            ind = dict_num(bra_basis, N_e);
                            auto it = hashT.find(ind);
                            if(it!=hashT.end())
                            {
                                indices_thread[count] = it->second;
                                data_thread[count] = A_list[m1*N_o*(2*N_o-1)+m2*(2*N_o-1)+m+N_o-1] * A_matirx[j1] * B_matirx[j2] * phase;
                                ++count;
                            }
                            #ifdef __DEBUG_D00A_GENERATE
                            printf("\thopping");
                            printf(" %d->%d ", m1+layerm1*N_o, m1+m+layerm1m*N_o);
                            printf("and %d->%d(m1=%d, m2=%d, m=%d):\t", m2+layerm2*N_o, m2-m+layerm2m*N_o, m1, m2,m);
                            print_basis(bra_basis, N_e, N_o);
                            printf("(j1=%d,j2=%d)(layerm1=%d,layerm1m=%d,layerm2=%d,layerm2m=%d)",j1,j2,layerm1,layerm1m,layerm2,layerm2m);
                            if(it!=hashT.end())
                                printf("\t resu num: %lu, ind: %lu, phase=%+d, resu=%lf\n", ind, it->second, phase, data_thread[count-1]);
                            else
                                printf("\t resu num: %lu, ind: hashT.end()\n", ind);
                            #endif
                        }
                    }
                }
            }
        }
    }

    free(bin_basis);
    free(bra_basis);
    return count;
}

extern "C" void D00AOp(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, int num_thread)
{
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        D00AOpGenerate_thread(A_matirx, B_matirx, A_list, vec, resu, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);
}

//#define __DEBUG_D00A_GENERATE

void D00AOpGenerate_thread(const Data_dtype *A_matirx, const Data_dtype *B_matirx, const Data_dtype *A_list, 
                            const Data_dtype *vec, Data_dtype *resu, Data_basis *basis, const hash & hashT, 
                            size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    int layers = 2;
    const Data_basis *basis_t = basis + row * N_e;
    size_t ind, N_orbit = layers * N_o;
    char phase, cphase;
    double me;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
    
    #ifdef __DEBUG_D00A_GENERATE
    printf("row=%lu\t", row);print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    resu[row] = 0;
    for(int m=1-(int)N_o;m<(int)N_o;++m)
    {
        for(int j1=0;j1<layers*layers;j1++)
        {
            for(int j2=0;j2<layers*layers;j2++)
            {
                if( std::abs(A_matirx[j1])>1e-15 && std::abs(B_matirx[j2])>1e-15 )
                {
                    int layerm1 = j1/2;
                    int layerm1m = j1&1;
                    int layerm2 = j2/2;
                    int layerm2m = j2&1;
                    for(int m1=0;m1<(int)N_o;++m1)
                    {
                        for(int m2=0;m2<(int)N_o;++m2)
                        {
                            me = A_list[m1*N_o*(2*N_o-1)+m2*(2*N_o-1)+m+N_o-1];
                            if( std::abs(me)<1e-15 )
                                continue;
                            basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
                            if(!hopping_LayertoLayer(bin_basis, m2, m, layerm2, layerm2m, N_e, N_o, cphase))
                                continue;


                            phase = cphase;
                            if(!hopping_LayertoLayer(bin_basis, m1, -m, layerm1, layerm1m, N_e, N_o, phase))
                                continue;
                            phase *= cphase;

                            binary_to_basis(bin_basis, bra_basis, N_orbit);

                            ind = dict_num(bra_basis, N_e);
                            auto it = hashT.find(ind);
                            if(it!=hashT.end())
                                resu[row] += vec[it->second] * me * A_matirx[j1] * B_matirx[j2] * phase;
                            #ifdef __DEBUG_D00A_GENERATE
                            printf("\thopping");
                            printf(" %d->%d ", m1+layerm1*N_o, m1+m+layerm1m*N_o);
                            printf("and %d->%d(m1=%d, m2=%d, m=%d):\t", m2+layerm2*N_o, m2-m+layerm2m*N_o, m1, m2,m);
                            print_basis(bra_basis, N_e, N_o);
                            printf("(j1=%d,j2=%d)(layerm1=%d,layerm1m=%d,layerm2=%d,layerm2m=%d)",j1,j2,layerm1,layerm1m,layerm2,layerm2m);
                            if(it!=hashT.end())
                                printf("\t resu num: %lu, ind: %lu, phase=%+d, resu=%lf\n", ind, it->second, phase, data_thread[count-1]);
                            else
                                printf("\t resu num: %lu, ind: hashT.end()\n", ind);
                            #endif
                        }
                    }
                }
            }
        }
    }

    free(bin_basis);
    free(bra_basis);
    return;
}


extern "C" void nl0AOp(Data_dtype *A_matirx, Data_dtype *vec, Data_dtype *resu, int l, int num_thread)
{
    struct CG2 cg;
    cg.cg = NULL;
    int st = 0;
    st = createCG2Bdoy(&cg, A.N_o);
    if(st!=0)
    {
        printf("Dl0AOp\t%d\n", st);
        exit(st);
    }
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        nl0AOpGenerate_thread(A_matirx, l, cg, vec, resu, A.basis, A.hashT, A.dim, A.N_o, A.N_e, i);
    destoryCG2(&cg);
}

//#define __DEBUG_NL0A_Op

void nl0AOpGenerate_thread(const Data_dtype *A_matirx, int l, const struct CG2 &cg, const Data_dtype *vec, Data_dtype *resu, Data_basis *basis, const hash & hashT, 
                            size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    int layers = 2;
    const Data_basis *basis_t = basis + row * N_e;
    size_t ind, N_orbit = layers * N_o;
    char phase;
    double me;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);
    
    #ifdef __DEBUG_NL0A_Op
    printf("row=%lu\t", row);print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    resu[row] = 0;
    for(int m1=0;m1<N_o;++m1)
    {
        for(int j=0;j<layers*layers;j++)
        {
            if( std::abs(A_matirx[j])>1e-15 )
            {
                int layerm1 = j/2;
                int layerm1m = j&1;
                if(!hopping_LayertoLayer(bin_basis, m1, 0, layerm1, layerm1m, N_e, N_o, phase))
                    continue;
                auto me = cg.cg[l][(N_o-1-m1)*N_o+(N_o-1-m1)] * cg.cg[l][0]; // m1 = s -> N_o-1
                me *= ((N_o-1+l+m1)&1)?-1:1;

                binary_to_basis(bin_basis, bra_basis, N_orbit);
                ind = dict_num(bra_basis, N_e);
                auto it = hashT.find(ind);
                if(it!=hashT.end())
                    resu[row] += vec[it->second] * me * A_matirx[j] * phase;

                #ifdef __DEBUG_NL0A_Op
                printf("\thopping");
                printf(" %d->%d(m1=%d)\t", m1+layerm1*N_o, m1+m+layerm1m*N_o);
                print_basis(bra_basis, N_e, N_o);
                printf("(j=%d)(layerm1=%d,layerm1m=%d)",j,layerm1,layerm1m);
                if(it!=hashT.end())
                    printf("\t resu num: %lu, ind: %lu, phase=%+d, resu=%lf\n", ind, it->second, phase, data_thread[count-1]);
                else
                    printf("\t resu num: %lu, ind: hashT.end()\n", ind);
                #endif
                restore_LayertoLayer(bin_basis, m1, 0, layerm1, layerm1m, N_o);
            }
        }
    }
    resu[row] *= (double)N_o/std::sqrt(4.0*M_PI*(2*l+1));

    free(bin_basis);
    free(bra_basis);
    return;
}
