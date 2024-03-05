

// FermionHamiltonianSymmetry.c

extern "C" size_t twoBodyCount(MKL_INT *indptr, int num_thread)
{
    indptr[0] = 0;

    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        indptr[i+1] = twoBodyCount_thread(A.basis, A.dim, A.N_o, A.N_e, i);

    for(size_t i=1;i<=A.dim;++i)
        indptr[i] += indptr[i-1];
    return indptr[A.dim];
}

size_t twoBodyCount_thread(Data_basis *k1_basis, size_t k1_dim, \
                    size_t N_o, size_t N_e, size_t i)
{
    size_t count = 0;
    Data_basis *basis_t = k1_basis + i * N_e;
    size_t j1, j2, j3, j4, J12, J34;
    size_t N_orbit = 2 * N_o;

    // CDW term
    //cdw_energy(k1_basis+i*N_e, N_phi, N_e);
    //printf("\t%u %lu\n", index[dict_num(basis_t, N_e)], i);
    count += tight_bond_count(basis_t, N_o, N_e);
    //printf("row=%lu, count=%lu\n", i, count);

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    Data_basis *diff = (Data_basis *)malloc((N_orbit-N_e) * sizeof(Data_basis));
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
            ++k;
        }
    }
    basis_t = k1_basis + i * N_e;

    for(size_t u=0;u<N_e-1;++u)
    {
        j3 = basis_t[u];
        for(size_t v=u+1;v<N_e;++v)
        {
            j4 = basis_t[v];
            J34 = j3 + j4;
            for(size_t s=0;s<N_orbit-N_e-1;++s)
            {
                j2 = diff[s];
                if(layerQ(j2, j3, N_o)==0)
                    continue;
                for(size_t t=s+1;t<N_orbit-N_e;++t)
                {
                    j1 = diff[t];
                    if(layerQ(j1, j4, N_o)==0)
                        continue;
                    J12 = j1 + j2;
                    if(J12!=J34)
                        continue;
                    ++count;
                    //printf("row=%lu, count=%lu\n", i, count);
                    //printf("\tFQHE:\t%d %d %d %d \ttype:\t%d\n", j1 ,j2, j3, j4, layerQ(j1, j2, (int)N_phi));
                }
            }
        }
    }
    free(bra_i);
    free(bra_j);
    free(diff);
    return count+1;
}

size_t tight_bond_count(Data_basis *basis, size_t N_phi, size_t N_e)
{
    size_t j1, j2, count = 0;
    Data_basis *temp = (Data_basis *)calloc(N_phi*2, sizeof(Data_basis));
    for(size_t i=0;i<N_e;++i)
        temp[basis[i]] = 1;

    for(size_t i=0;i<N_e;++i)
    {
        j2 = basis[i];
        j1 = tight_up(j2, N_phi);
        if(temp[j1]==0)
        {
            ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
        }
        /*j1 = tight_down(j2, N_phi);
        if(temp[j1]==0)
        {
            ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
        }*/
    }

    free(temp);
    return count;
}

size_t tight_up(size_t j1, size_t N_phi)
{
    if(j1>=N_phi)
        return j1 - N_phi;
    else
        return j1 + N_phi;
}

/*
size_t tight_down(size_t j1, size_t N_phi)
{
    if(j1<N_phi)
        return j1 + N_phi;
    else
        return j1 - N_phi;
}*/

size_t layerQ(size_t j1, size_t j2, size_t N_phi)
{
    return (j1<N_phi)?((j2<N_phi)?1:0):((j2<N_phi)?0:1);
}


extern "C" size_t twoBodyGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double *cdw, double t, int num_thread)
{
    //for(size_t i=0;i<2*N_phi-1;++i)
    //{
    //  for(size_t j=0;j<2*N_phi-1;++j)
    //  {
    //      printf("%lf+I%lf\t", creal(A1_list[i*(2*N_phi-1)+j]), cimag(A1_list[i*(2*N_phi-1)+j]));
    //  }printf("\n");
    //}
    indptr[0] = 0;
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
    {
        auto count = twoBodyGenerate_thread(indptr, indices, data, A1_list, A12_list, cdw, \
                                    A.basis, A.hashT, t, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function twoBodyGenerate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

size_t twoBodyGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double *cdw, \
                    Data_basis *k1_basis, const hash & hashT, double t, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    size_t bra_num, ind, count = 0;
    Data_basis *basis_t = k1_basis + row * N_e;
    size_t j1, j2, j3, j4, j13, j14, j23, j24, J12, J34, J1, J2, J3, J4;
    size_t N_orbit = 2 * N_o;
    Data_dtype dig = 0.0;
    MKL_INT *indices_thread = indices + indptr[row];
    Data_dtype* data_thread = data + indptr[row];
    int phase;

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    Data_basis *diff = (Data_basis *)malloc((N_orbit-N_e) * sizeof(Data_basis));
    // CDW term
    //cdw_energy(k1_basis+i*N_e, N_phi, N_e); //dig
    //printf("\t%u %lu\n", index[dict_num(basis_t, N_e)], i);
    
    // tight bond term
    count = tight_bond(basis_t, bra_i, bra_j, hashT, indices_thread, data_thread, \
                        t, N_o, N_e, count);// update count
    //printf("row=%lu, count=%lu\n", row, count);
    // FQHE term and CDW term
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
    basis_t = k1_basis + row * N_e;

    for(size_t u=0;u<N_e-1;++u)
    {
        j3 = basis_t[u];
        for(size_t v=u+1;v<N_e;++v)
        {
            j4 = basis_t[v];
            J34 = j3 + j4;
            for(size_t s=0;s<N_orbit-N_e-1;++s)
            {
                j2 = diff[s];
                if(layerQ(j2, j3, N_o)==0)
                    continue;
                for(size_t t=s+1;t<N_orbit-N_e;++t)
                {
                    j1 = diff[t];
                    //if(row==0)
                    //  printf("j!%d, %d, %d, %d\n", j1, j2, j3, j4);
                    if(layerQ(j1, j4, N_o)==0)
                        continue;
                    J12 = j1 + j2;
                    if(J12!=J34)
                        continue;
                    bra_num = intra_hopping(basis_t, bra_i, bra_j, 2*N_o, N_e, j1, j2, j3, j4, &phase);
                    auto it = hashT.find( bra_num );
                    if(it==hashT.end()) 
                    {   
                        printf("%lu %lu %lu %lu %lu %lu\n", j1, j2, j3, j4, bra_num, row);
                        printf("twoBodyGenerate_thread: Not found\n");
                        continue;
                    }
                    ind = it->second;
                    //ind = index[bra_num];
                    J1 = j1%N_o;
                    J2 = j2%N_o;
                    J3 = j3%N_o;
                    J4 = j4%N_o;
                    if(layerQ(j1, j2, N_o)==0)
                        dig = 2.0 * A12_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4] * phase;
                    else
                        dig = 2.0 * (A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4]-A1_list[J2*N_o*N_o*N_o + J1*N_o*N_o + J3*N_o + J4]) * phase;
                    //printf("%d, %d, %d, %d, %d, %d, %d, %d %lu %lu %d\n", \
                                j1, j2, j3, j4, J1, J2, J3, J4, ind, row, layerQ(j1, j2, N_o));
                    //printf("%lf+I%lf, phase=%d\n", real(dig), imag(dig), phase);
                    //printf("%lf, %lf, %lf, %lu\n", A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4], A1_list[J2*N_o*N_o*N_o + J1*N_o*N_o + J3*N_o + J4], A1_list[17], u);

                    /*if(layerQ(j1, j2, N_o)==0)    //inter-layer
                    {
                        bra_num = inter_hopping(basis_t, bra_i, bra_j, 2*N_o, N_e, j1, j2, j3, j4, &phase);
                        ind = index[bra_num];
                        J1 = j1;
                        J2 = j2;
                        J3 = j3;
                        J4 = j4;
                        while(J1>=N_o)
                            J1 -= N_o;
                        while(J2>=N_o)
                            J2 -= N_o;
                        while(J3>=N_o)
                            J3 -= N_o;
                        while(J4>=N_o)
                            J4 -= N_o;
                        //j13 = N_o-1+J1-J3;
                        //j14 = N_o-1+J1-J4;
                        //j23 = j13-J1+J2;
                        //j24 = j14-J1+J2;
                        //dig = A12_list[j13*(2*N_o-1)+j14] * phase;
                        dig = A12_list[J1*u*u*u + J2*u*u + J3*u + J4] * phase;
                    }
                    else                                //intra-layer
                    {
                        j13 = N_phi-1+j1-j3;
                        j14 = j13+j3-j4;
                        j23 = j13-j1+j2;
                        j24 = j14-j1+j2;
                        bra_num = intra_hopping(basis_t, bra_i, bra_j, 3*N_phi, N_e, j1, j2, j3, j4, &phase);
                        ind = index[bra_num];
                        dig = 2.0 * (A1_list[j13*(2*N_phi-1)+j14]-A1_list[j23*(2*N_phi-1)+j24]) * phase;
                    
                        //printf("%d %d %d %d %lu %lu %lu\n", \
                                j1, j2, j3, j4, ind, bra_num, row);
                        //printf("%lf+I%lf\n", creal(dig), cimag(dig));
                    }*/
                    indices_thread[count] = ind;
                    data_thread[count] = dig;
                    ++count;
                    //printf("row=%lu, count=%lu\n", row, count);
                }
            }
        }
    }
    dig *= 0.0;
    for(size_t s1=0;s1<N_e;++s1)
    {
        for(size_t s2=s1+1;s2<N_e;++s2)
        {
            j3 = basis_t[s1];
            j4 = basis_t[s2];
            J3 = j3%N_o;
            J4 = j4%N_o;
            j1 = j4;
            j2 = j3;
            J1 = J4;
            J2 = J3;
            if(layerQ(j1, j2, N_o)==0)// inter
                dig += 2.0 * A12_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4];
            else// intra
                dig += 2.0 * (A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4]-A1_list[J2*N_o*N_o*N_o + J1*N_o*N_o + J3*N_o + J4]);
            //printf("%d, %d, %d, %d, %lu\n", j1, j2, j3, j4, row);
            //printf("%lf+I%lf\n", real(dig), imag(dig));
            //printf("%d %d %d\n",layerQ(j3, j4, (int)N_phi),((j2<2*N_phi) && (j2>=N_phi)), j2>=N_phi);
            /*if(layerQ(j3, j4, (int)N_phi)==0) //inter-layer
            {
                J1 = j1;
                J2 = j2;
                J3 = j3;
                J4 = j4;
                while(J1>=N_phi)
                    J1 -= N_phi;
                while(J2>=N_phi)
                    J2 -= N_phi;
                while(J3>=N_phi)
                    J3 -= N_phi;
                while(J4>=N_phi)
                    J4 -= N_phi;
                j13 = N_phi-1+J1-J3;
                j14 = N_phi-1+J1-J4;
                j23 = j13-J1+J2;
                j24 = j14-J1+J2;
                //printf("%lf+I%lf\n", creal(dig), cimag(dig));
                dig += A12_list[j13*(2*N_phi-1)+j14];
                //printf("%d, %d, %d, %d, %d, %d, %d, %d %lu %lu %lu\n", \
                                    j1, j2, j3, j4, j13, j14, j23, j24, ind, bra_num, row);
                //printf("%lf+I%lf\n", creal(dig), cimag(dig));
            }
            else if(s2>s1)                          //intra-layer
            {
                j13 = N_phi-1+j1-j3;
                j14 = j13+j3-j4;
                j23 = j13-j1+j2;
                j24 = j14-j1+j2;
                //printf("%lf+I%lf\n", creal(dig), cimag(dig));
                dig += 2 * (A1_list[j13*(2*N_phi-1)+j14]-A1_list[j23*(2*N_phi-1)+j24]);
                //printf("%d, %d, %d, %d, %d, %d, %d, %d %lu %lu %lu\n", \
                                    j1, j2, j3, j4, j13, j14, j23, j24, ind, bra_num, row);
                //printf("%lf+I%lf\n", creal(dig), cimag(dig));
            }*/
        }
    }
            //printf("%d, %d, %d, %d, %lu\n", j1, j2, j3, j4, row);
    //printf("%lu %lu %u %lf+I%lf\n", row, count, indptr[row+1], real(dig), imag(dig));
    //getchar();
    indices_thread[count] = row;
    data_thread[count] = dig + cdw_energy(basis_t, cdw, N_o, N_e);
    //printf("%lu %lf %lf+I%lf\n", row, cdw_energy(basis_t, cdw, N_phi, N_e), creal(dig), cimag(dig));
    ++count;


    free(bra_i);
    free(bra_j);
    free(diff);
    return count;
}

double cdw_energy(Data_basis *basis, double *cdw, size_t N_phi, size_t N_e)
{
    int level[2] = {0, 0};
    for(size_t i=0;i<N_e;++i)
    {
        //printf("%d ", basis[i]);
        if(basis[i]<N_phi)
            ++level[0];
        else //if(basis[i]<2*N_phi)
            ++level[1];
    }
    //printf("\t%d %d %d:", level[0], level[1], level[2]);
    return cdw[0]*level[0] + cdw[1]*level[1];
}

size_t tight_bond(Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
                const hash & hashT, MKL_INT *indices_thread, Data_dtype *data_thread, \
                double t, size_t N_phi, size_t N_e, size_t count)
{
    size_t j1, j2, ind;
    int phase;
    Data_basis *temp = (Data_basis *)calloc(N_phi*2, sizeof(unsigned));
    for(size_t i=0;i<N_e;++i)
        temp[basis_t[i]] = 1;

    for(size_t i=0;i<N_e;++i)
    {
        j2 = basis_t[i];
        j1 = tight_up(j2, N_phi);
        if(temp[j1]==0)
        {
            ind = tight_bond_hopping(basis_t, bra_i, bra_j, N_e, N_phi, j1, j2, &phase);
            //indices_thread[count] = index[ind];
            auto it = hashT.find( ind );
            if(it!=hashT.end())
            {
                indices_thread[count] = it->second;
                data_thread[count] = t * phase;
                ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
            }
            else
            {
                printf("tight_bond: Not found!\n");
            }
        }
        /*j1 = tight_down(j2, N_phi);
        if(temp[j1]==0)
        {
            ind = tight_bond_hopping(basis_t, bra_i, bra_j, N_e, N_phi, j1, j2, &phase);
            indices_thread[count] = index[ind];
            data_thread[count] = t * phase;
            ++count;//printf("\ttight_bond:\t%lu %lu\n", j1 ,j2);
        }*/
    }

    free(temp);
    return count;
}

size_t tight_bond_hopping(Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
                            size_t N_e, size_t N_phi, size_t j1, size_t j2, int *phase)
{
    size_t N_orbit = 2*N_phi;

    // initialization
    memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
    memset(bra_j, 0, (size_t)N_e * sizeof(Data_basis));

    // translation occupy -> binary
    for(size_t i=0;i<N_e;++i)
        ++bra_i[basis_t[i]];
    *phase = 0;

    size_t b = (j1>j2)?j2:j1;
    size_t e = (j1>j2)?j1:j2;

    //phase
    for(size_t i=b+1;i<e;++i)
        *phase += bra_i[i];
    // annihilation
                //for(auto ii=0;ii<N_orbit;++ii)
                //  printf("%u ", bra_i[ii]);
                //printf("\n");
    --bra_i[j2];
    // creation
    ++bra_i[j1];
                //for(auto ii=0;ii<N_orbit;++ii)
                //  printf("%u ", bra_i[ii]);
                //printf("\n");

    *phase = (*phase&1)?-1:1;
    
    // translation binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(bra_i[i])
        {
            bra_j[k] = i;
            ++k;
        }
    }

    return dict_num(bra_j, N_e);
}
