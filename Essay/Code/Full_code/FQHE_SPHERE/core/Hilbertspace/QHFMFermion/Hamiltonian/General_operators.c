

// daggerQ: 
//          0: a_{m-m_1}a_{m_1}  
//          1: a_{m_1}^daggera_{m-m_1}^dagger
extern "C" void pairingOp(Data_dtype *A_list, Data_dtype *A_matirx, Data_dtype *vec, Data_dtype *resu, int l, int m, 
                            int daggerQ, int num_thread)
{
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<B.dim;++i)
        pairingOpGenerate_thread(A_list, A_matirx, vec, resu, m, A.basis, A.hashT, A.dim, A.N_o, A.N_e,
                         B.basis, B.hashT, B.dim, B.N_o, B.N_e, daggerQ, i);
    //destoryCG2(&cg);
}

#ifndef __DEBUG_PAIRINGOPGENERATE
//#define __DEBUG_PAIRINGOPGENERATE
#endif

//a_{m-m_1}a_{m_1} 
void pairingOpGenerate_thread(const Data_dtype *A_list, const Data_dtype *A_matirx, const Data_dtype *vec, Data_dtype *resu, int m,
                            Data_basis *rbasis, const hash & rhashT, size_t rk1_dim, size_t rN_o, size_t rN_e, 
                            Data_basis *lbasis, const hash & lhashT, size_t lk1_dim, size_t lN_o, size_t lN_e, int daggerQ, size_t row)
{
    int layers = 2;
    const Data_basis *lbasis_t = lbasis + row * lN_e;
    size_t N_o = lN_o;
    size_t N_orbit = layers * N_o;
    long int ind;
    int phase;
    double me;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(lN_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(lbasis_t, bin_basis, N_orbit, lN_e);

    Data_basis *ket = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra = (Data_basis *)malloc(rN_e * sizeof(Data_basis));
    //Data_basis *diff = (Data_basis *)malloc((N_orbit-rN_e) * sizeof(Data_basis));

    #ifdef __DEBUG_PAIRINGOPGENERATE
    printf("row=%lu\t", row);print_basis(lbasis_t, lN_e, N_o);printf("\n");
    #endif

    int temp, j1, j2;
    double comm;

    resu[row] = 0;
    for(int m1=0;m1<N_o;++m1)
    {
        j1 = m+N_o-1-m1;
        if(j1<0||j1>=N_o)
            continue;
        for(int j=0;j<layers*layers;j++)
        {
            if( std::abs(A_matirx[j])>1e-15 )
            {
                int layerm1 =  j/2; //A^dagger
                int layerm1m = j&1;
                j1 = m+N_o-1-m1;
                j2 = m1;

                comm = A_list[j2*N_o+j1];
                if(daggerQ==0)
                {
                    std::swap(layerm1, layerm1m); // A^dagger
                    std::swap(j1, j2);
                    ind = creation_two(bin_basis, bra, ket, N_o, lN_e, j1+layerm1*N_o, j2+layerm1m*N_o, phase);
                }
                else
                {
                    ind = destroy_two(bin_basis, bra, ket, N_o, lN_e, j1+layerm1*N_o, j2+layerm1m*N_o, phase);
                }

                if(ind<0)
                    continue;
                auto it = rhashT.find(ind);
                if(it!=rhashT.end())
                    resu[row] += vec[it->second] * comm * A_matirx[j] * phase;//A_matirx[layerm1*layers+layerm1m] * phase;

                #ifdef __DEBUG_PAIRINGOPGENERATE
                printf("\thopping");
                if(daggerQ==0)
                    printf(" c^†_%d c^†_%d (j1=%d, j2=%d, m1=%d, m=%d)\t", j1+layerm1*N_o, j2+layerm1m*N_o, j1, j2, m1, m);
                else
                    printf(" c_%d c_%d (j1=%d, j2=%d, m1=%d, m=%d)\t", j1+layerm1*N_o, j2+layerm1m*N_o, j1, j2, m1, m);
                print_basis(bra, rN_e, N_o);
                printf("(j=%d)(layerm1=%d,layerm1m=%d)",j,layerm1,layerm1m);
                if(it!=rhashT.end())
                    printf("\t resu num: %ld, ind: %lu, phase=%+d, A_list=%lf\n", ind, it->second, phase, A_list[(m+N_o-1-m1)*N_o+m1]);
                else
                    printf("\t resu num: %ld, ind: hashT.end()\n", ind);
                //getchar();
                #endif
            }
        }
    }
    //getchar();

    free(bra);
    free(ket);
    //free(diff);
    free(bin_basis);
    free(bra_basis);
    return;
}

#ifndef __DEBUG_CREATION_TWO
//#define __DEBUG_CREATION_TWO
#endif

// a_{j1}^dagger a_{j2}^dagger |>
long creation_two(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e, int j1, int j2, int& phase)
{
    size_t N_orbit = 2*N_o;
    memset(ket, 0, (size_t)N_orbit * sizeof(Data_basis));
    if(bin_basis[j1]||bin_basis[j2]||j1==j2)
        return -1;
    for(int i=0;i<N_orbit;++i)
        ket[i] = bin_basis[i];

    #ifdef __DEBUG_CREATION_TWO
    printf("\t");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    #endif

    ket[j1] = 1;
    ket[j2] = 1;

    #ifdef __DEBUG_CREATION_TWO
    printf(" --> ");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    printf("\n");
    #endif

    phase = 0;
    for(int i=std::min(j1,j2)+1;i<std::max(j1,j2);++i)
        phase += ket[i];
    phase = (phase&1)?-1:1;
    if(j2<j1)
        phase *= -1;

    // translation binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(ket[i])
        {
            bra[k] = i;
            ++k;
        }
    }

    return dict_num(bra, N_e+2);
}

// a_{j1} a_{j2} |>
long destroy_two(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e, int j1, int j2, int& phase)
{
    size_t N_orbit = 2*N_o;
    memset(ket, 0, (size_t)N_orbit * sizeof(Data_basis));
    if(bin_basis[j1]==0||bin_basis[j2]==0||j1==j2)
        return -1;
    for(int i=0;i<N_orbit;++i)
        ket[i] = bin_basis[i];

    #ifdef __DEBUG_CREATION_TWO
    printf("\t");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    #endif

    ket[j1] = 0;
    ket[j2] = 0;

    #ifdef __DEBUG_CREATION_TWO
    printf(" --> ");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    printf("\n");
    #endif

    phase = 0;
    for(int i=std::min(j1,j2)+1;i<std::max(j1,j2);++i)
        phase += ket[i];
    phase = (phase&1)?-1:1;
    if(j2>j1)
        phase *= -1;

    // translation binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(ket[i])
        {
            bra[k] = i;
            ++k;
        }
    }

    return dict_num(bra, N_e-2);
}


extern "C" void pairingABlmOp(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, 
                                int l, int m, int la, int lb, int num_thread)
{
    if(m!=0)
    {
        Data_basis * lbasis = B.basis;
        const hash & lhashT = B.hashT;
        size_t lk1_dim = B.dim;
        # pragma omp parallel for num_threads(num_thread)
        for(size_t i=0;i<B.dim;++i)
            pairingABlmOpGenerate_thread(A_list, A_matirx, B_matirx, vec, resu, l, m, la, lb, 
                                        A.basis, A.hashT, A.dim, A.N_o, A.N_e, 
                                        lbasis, lhashT, lk1_dim, i);
    }
    else
    {
        Data_basis * lbasis = A.basis;
        const hash & lhashT = A.hashT;
        size_t lk1_dim = A.dim;
        # pragma omp parallel for num_threads(num_thread)
        for(size_t i=0;i<lk1_dim;++i)
            pairingABlmOpGenerate_thread(A_list, A_matirx, B_matirx, vec, resu, l, m, la, lb, 
                                        A.basis, A.hashT, A.dim, A.N_o, A.N_e, 
                                        lbasis, lhashT, lk1_dim, i);
    }
}

#ifndef __DEBUG_PAIRINGABLMOP
//#define __DEBUG_PAIRINGABLMOP
#endif

void pairingABlmOpGenerate_thread(const Data_dtype *A_list, const Data_dtype *A_matirx,  const Data_dtype *B_matirx, 
                            const Data_dtype *vec, Data_dtype *resu, int l, int m, int la, int lb,
                            Data_basis *rbasis, const hash & rhashT, size_t rk1_dim, size_t N_o, size_t N_e, 
                            Data_basis *lbasis, const hash & lhashT, size_t lk1_dim, size_t row)
{
    int layers = 2;
    const Data_basis *basis_t = lbasis + row * N_e;
    size_t ind, N_orbit = layers * N_o;
    int phase;
    double dig, me, Aab, Bcd;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);

    Data_basis *ket = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    Data_basis *diff = (Data_basis *)malloc((N_orbit-N_e+2) * sizeof(Data_basis));
    
    #ifdef __DEBUG_PAIRINGABLMOP
    printf("row=%lu\t", row);print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    size_t k = 0;
    for(size_t j=0,ll=0;j<N_orbit;++j)
    {
        if((*basis_t == j)&&(ll<N_e))
        {
            ++basis_t;
            ++ll;
        }
        else
        {
            diff[k] = j;
            //printf("%lu ", j);
            ++k;
        }
    }//printf("\n");
    basis_t = lbasis + row * N_e;

    int j1, j2, j3, j4, J34, J12, J1, J2, J3, J4;

    for(size_t u=0;u<N_e;++u)
    {
        j1 = basis_t[u];
        J1 = j1%N_o;
        for(size_t v=0;v<N_e;++v)
        {
            if(u==v) continue;
            j2 = basis_t[v];
            J2 = j2%N_o;
            J12 = J1 + J2;
            if(std::abs(J12-int(N_o)+1)>la+0.1 ) continue;
            Aab = checkLayer(A_matirx, j2, j1, N_o, 2, false);
            if(std::abs(Aab)<1e-15) continue;
            diff[N_orbit-N_e] = j1;
            diff[N_orbit-N_e+1] = j2;
            for(size_t s=0;s<N_orbit-N_e+2;++s)
            {
                j3 = diff[s];
                J3 = j3%N_o;
                for(size_t t=0;t<N_orbit-N_e+2;++t)
                {
                    if(s==t) continue;
                    j4 = diff[t];
                    J4 = j4%N_o;
                    J34 = J3 + J4;
                    if(std::abs(J34-int(N_o)+1)>lb+0.1 ) continue;
                    if(J12!=J34+m)
                        continue;
                    Bcd = checkLayer(B_matirx, j4, j3, N_o, 2, true);// B^\dagger
                    if(std::abs(Bcd)<1e-15) continue;

                    ind = pairing_hopping(bin_basis, bra, ket, N_o, N_e, j1, j2, j3, j4, phase);
                    if(ind<0) continue;
                    auto it = rhashT.find(ind);
                    if(it!=rhashT.end())
                        resu[row] += vec[it->second] * Aab * Bcd * phase * A_list[J1*N_o*N_o*N_o+J2*N_o*N_o+J3*N_o+J4];

                    #ifdef __DEBUG_PAIRINGABLMOP
                    printf("\thopping");
                    printf(" c^†_%d c^†_%d c_%d c_%d (j4=%d, j3=%d, j2=%d, j1=%d)\t", j4, j3, j2, j1, j4, j3, j2, j1);
                    print_basis(bra, N_e, N_o);
                    if(it!=rhashT.end())
                        printf("\t resu num: %ld, ind: %lu(%lf), phase=%+d, A_list=%lf, Aab=%lf, Bcd=%lf: %lf\n", 
                            ind, it->second, vec[it->second], phase, A_list[J1*N_o*N_o*N_o+J2*N_o*N_o+J3*N_o+J4], Aab, Bcd,resu[row]);
                    else
                        printf("\t resu num: %ld, ind: hashT.end()\n", ind);
                    #endif
                }
            }
        }
    }
    // diagonal term
    /*if(m==0)
    {
        for(size_t u=0;u<N_e;++u)
        {
            j1 = basis_t[u]; 
            J1 = j1%N_o;
            for(size_t v=0;v<N_e;++v)
            {
                if(u==v) continue;
                j2 = basis_t[v];
                J2 = j2%N_o;
                J12 = J1 + J2;
                printf("j1=%d, j2=%d ", j1, j2);
                if(std::abs(J12-int(N_o)+1)>la+0.1 ) 
                {
                    printf("J12=%d, std::abs(J12-int(N_o)+1)=%d, la=%d\n", J12, std::abs(J12-int(N_o)+1), la);
                    continue;
                }
                Aab = checkLayer(A_matirx, j2, j1, N_o, 2, false);
                if(std::abs(Aab)<1e-15)
                {
                    printf("Aab=%lf\n", Aab);
                    continue;
                }
                j3 = j1; j4 = j2;
                J3 = J1; J4 = J2;
                for(int s=0;s<2;s++)
                {
                    Bcd = checkLayer(B_matirx, j4, j3, N_o, 2, true);
                    if(std::abs(Bcd)<1e-15)
                    {
                        printf("Bcd=%lf\n", Aab);
                        continue;
                    }
                    ind = pairing_hopping(bin_basis, bra, ket, N_o, N_e, j1, j2, j3, j4, phase);
                    resu[row] += vec[row] * Aab * Bcd * phase * A_list[J1*N_o*N_o*N_o+J2*N_o*N_o+J3*N_o+J4];

                    #ifdef __DEBUG_PAIRINGABLMOP
                    printf("\tdiag");
                    printf(" c^†_%d c^†_%d c_%d c_%d (j4=%d, j3=%d, j2=%d, j1=%d)\t", j4, j3, j2, j1, j4, j3, j2, j1);
                    print_basis(bra, N_e, N_o);
                    auto it = rhashT.find(ind);
                    if(it!=rhashT.end())
                        printf("\t resu num: %ld, ind: %lu(%lf), phase=%+d, A_list=%lf, Aab=%lf, Bcd=%lf: %lf\n", 
                            ind, it->second, vec[it->second], phase, A_list[J1*N_o*N_o*N_o+J2*N_o*N_o+J3*N_o+J4], Aab, Bcd,resu[row]);
                    else
                        printf("\t resu num: %ld, ind: hashT.end()\n", ind);
                    #endif

                    std::swap(j3, j4); std::swap(J3, J4);
                }
            }
        }
    }*/
    #ifdef __DEBUG_PAIRINGABLMOP
    printf("\t Finally: %lf\n", resu[row]);
    getchar();
    #endif

    free(bin_basis);
    free(bra_basis);
    free(bra);
    free(ket);
    free(diff);
    return;
}

Data_dtype checkLayer(const Data_dtype *A_matirx, int j1, int j2, int N_o, int layers, bool daggerQ)
{
    int alpha, beta;
    alpha = int(j1/N_o);
    beta  = int(j2/N_o);
    if(daggerQ) std::swap(alpha, beta);
    return A_matirx[alpha*layers+beta];
}

#ifndef __DEBUG_PAIRING_HOPPING
//#define __DEBUG_PAIRING_HOPPING
#endif

// c_{j4}^dagger c_{j3}^dagger c_{j2} c_{j1} |>
long int pairing_hopping(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e,
                         int j1, int j2, int j3, int j4, int& phase)
{
    size_t N_orbit = 2*N_o;
    memset(ket, 0, (size_t)N_orbit * sizeof(Data_basis));
    if(bin_basis[j1]==0||bin_basis[j2]==0||j1==j2)
        return -1;
    for(int i=0;i<N_orbit;++i)
        ket[i] = bin_basis[i];

    #ifdef __DEBUG_PAIRING_HOPPING
    printf("\t");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    #endif

    // destroy
    ket[j1] = 0;
    ket[j2] = 0;
    phase = 0;
    for(int i=std::min(j1,j2)+1;i<std::max(j1,j2);++i)
        phase += ket[i];
    if(j2<j1)
        phase += 1;

    #ifdef __DEBUG_PAIRING_HOPPING
    printf(" --> ");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    #endif

    // creation
    if(ket[j3]||ket[j4]||j3==j4)
        return -1;
    ket[j3] = 1;
    ket[j4] = 1;
    for(int i=std::min(j3,j4)+1;i<std::max(j3,j4);++i)
        phase += ket[i];
    if(j3<j4)
        phase += 1;
    phase = (phase&1)?-1:1;

    #ifdef __DEBUG_PAIRING_HOPPING
    printf(" --> ");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    printf("\n");
    #endif

    // translation binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(ket[i])
        {
            bra[k] = i;
            ++k;
        }
    }

    return dict_num(bra, N_e);
}



extern "C" void twoOpGeneral(Data_dtype *A_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, int m, int num_thread)
{
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<B.dim;++i)
        twoOpGeneralGenerate_thread(A_list, A_matirx, vec, resu, m,
                                    A.basis, A.hashT, A.dim, A.N_o, A.N_e, 
                                    B.basis, B.hashT, B.dim, i);
}

#ifndef __DEBUG_TWOOPGENERAL
//#define __DEBUG_TWOOPGENERAL
#endif

void twoOpGeneralGenerate_thread(const Data_dtype *A_list, const Data_dtype *A_matirx, const Data_dtype *vec, Data_dtype *resu, int m,
                            Data_basis *rbasis, const hash & rhashT, size_t rk1_dim, size_t N_o, size_t N_e, 
                            Data_basis *lbasis, const hash & lhashT, size_t lk1_dim, size_t row)
{
    int layers = 2;
    const Data_basis *basis_t = lbasis + row * N_e;
    size_t ind, N_orbit = layers * N_o;
    int phase;
    double dig, me, Aab;

    Data_basis *bin_basis = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_basis = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    memset(bin_basis, 0, (size_t)N_orbit * sizeof(Data_basis));
    basis_to_binary(basis_t, bin_basis, N_orbit, N_e);

    Data_basis *ket = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    Data_basis *diff = (Data_basis *)malloc((N_orbit-N_e+1) * sizeof(Data_basis));
    
    #ifdef __DEBUG_TWOOPGENERAL
    printf("row=%lu\t", row);print_basis(basis_t, N_e, N_o);printf("\n");
    #endif

    size_t k = 0;
    for(size_t j=0,ll=0;j<N_orbit;++j)
    {
        if((*basis_t == j)&&(ll<N_e))
        {
            ++basis_t;
            ++ll;
        }
        else
        {
            diff[k] = j;
            //printf("%lu ", j);
            ++k;
        }
    }//printf("\n");
    basis_t = lbasis + row * N_e;

    int j1, j2, J1, J2;

    for(size_t u=0;u<N_e;++u)
    {
        j1 = basis_t[u];
        J1 = j1%N_o;
        diff[N_orbit-N_e] = j1;
        for(size_t s=0;s<N_orbit-N_e+1;++s)
        {
            j2 = diff[s];
            J2 = j2%N_o;
            if(J1!=J2+m)
                continue;
            Aab = checkLayer(A_matirx, j2, j1, N_o, 2, true);// A^\dagger
            if(std::abs(Aab)<1e-15) continue;

            ind = nlm_hopping(bin_basis, bra, ket, N_o, N_e, j1, j2, phase);
            if(ind<0) continue;
            auto it = rhashT.find(ind);
            if(it!=rhashT.end())
                resu[row] += vec[it->second] * Aab * phase * A_list[J1*N_o+J2];

            #ifdef __DEBUG_TWOOPGENERAL
            printf("\thopping");
            printf(" c^†_%d c_%d (j2=%d, j1=%d)\t", j2, j1, j2, j1);
            print_basis(bra, N_e, N_o);
            if(it!=rhashT.end())
                printf("\t resu num: %ld, ind: %lu(%lf), phase=%+d, A_list=%lf, Aab=%lf: %lf\n", 
                    ind, it->second, vec[it->second], phase, A_list[J1*N_o+J2], Aab, resu[row]);
            else
                printf("\t resu num: %ld, ind: hashT.end()\n", ind);
            #endif
    }
    }
    #ifdef __DEBUG_TWOOPGENERAL
    printf("\t Finally: %lf\n", resu[row]);
    getchar();
    #endif

    free(bin_basis);
    free(bra_basis);
    free(bra);
    free(ket);
    free(diff);
    return;
}


#ifndef __DEBUG_NLM_HOPPING
//#define __DEBUG_NLM_HOPPING
#endif

// c_{j2}^dagger c_{j1} |>
long int nlm_hopping(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e,
                         int j1, int j2, int& phase)
{
    size_t N_orbit = 2*N_o;
    memset(ket, 0, (size_t)N_orbit * sizeof(Data_basis));
    if(bin_basis[j1]==0)
        return -1;
    for(int i=0;i<N_orbit;++i)
        ket[i] = bin_basis[i];

    #ifdef __DEBUG_NLM_HOPPING
    printf("\t");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    #endif

    // destroy
    ket[j1] = 0;

    #ifdef __DEBUG_NLM_HOPPING
    printf(" --> ");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    #endif

    // creation
    if(ket[j2])
        return -1;
    ket[j2] = 1;

    phase = 0;
    for(int i=std::min(j1,j2)+1;i<std::max(j1,j2);++i)
        phase += ket[i];
    phase = (phase&1)?-1:1;

    #ifdef __DEBUG_NLM_HOPPING
    printf(" --> ");
    for(int i=0;i<N_orbit;i++)
        printf("%d ", ket[i]);
    printf("\n");
    #endif

    // translation binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(ket[i])
        {
            bra[k] = i;
            ++k;
        }
    }

    return dict_num(bra, N_e);
}