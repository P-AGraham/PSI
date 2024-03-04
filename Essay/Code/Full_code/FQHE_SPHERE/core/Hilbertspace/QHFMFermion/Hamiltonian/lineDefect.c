

extern "C" size_t defectLineGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, double *A1_list, double *A12_list, double *cdw, double* defect_list, double t, int num_thread)
{
    indptr[0] = 0;
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
    {
        auto count = defectLineGenerate_thread(indptr, indices, data, A1_list, A12_list, cdw, defect_list, A.basis, A.hashT, t, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function defectLineGenerate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

size_t defectLineGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double *cdw, double* defect_list,\
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
            ++k;
        }
    }
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
                    indices_thread[count] = ind;
                    data_thread[count] = dig;
                    ++count;
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
        }
    }
    for(int i=0;i<N_o;++i)
        dig += defect_list[basis_t[i]%N_o]*((basis_t[i]<N_o)?1.0:-1.0);

    indices_thread[count] = row;
    data_thread[count] = dig + cdw_energy(basis_t, cdw, N_o, N_e);
    ++count;


    free(bra_i);
    free(bra_j);
    free(diff);
    return count;
}
