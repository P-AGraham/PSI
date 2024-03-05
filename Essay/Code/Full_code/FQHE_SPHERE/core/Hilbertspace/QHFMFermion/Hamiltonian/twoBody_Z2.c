

// FermionHamiltonianSymmetry.c

extern "C" size_t twoBodyZ2Count(MKL_INT *indptr, int num_thread)
{
    indptr[0] = 0;
    //indptr[A.dim] = 3300000000;
    //return 3300000000;

    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<A.dim;++i)
        indptr[i+1] = twoBodyZ2Count_thread(A.basis, A.hashT, A.Z2_phase, A.Z_2, A.dim, A.N_o, A.N_e, i);

    for(size_t i=1;i<=A.dim;++i)
        indptr[i] += indptr[i-1];
    return indptr[A.dim];
}

size_t twoBodyZ2Count_thread(Data_basis *k1_basis, const hash & hashT, char* Z2_phase, int Z_2, size_t k1_dim, \
                    size_t N_o, size_t N_e, size_t i)
{
    size_t count = 0;
    Data_basis *basis_t = k1_basis + i * N_e;
    size_t j1, j2, j3, j4, J12, J34, idx;
    size_t N_orbit = 2 * N_o;
    int phase, n_find;

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    Data_basis *diff = (Data_basis *)malloc((N_orbit-N_e) * sizeof(Data_basis));

    // CDW term
    //cdw_energy(k1_basis+i*N_e, N_phi, N_e);
    //printf("\t%u %lu\n", index[dict_num(basis_t, N_e)], i);
    count += tight_bond_Z2_count(basis_t, bra_i, bra_j, hashT, N_o, N_e);
    //printf("row=%lu, count=%lu\n", i, count);

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
                    intra_hopping(basis_t, bra_i, bra_j, 2*N_o, N_e, j1, j2, j3, j4, &phase);
                    idx = find_Z2_index(bra_j, hashT, N_e, N_o, n_find);
                    if(idx!=ULIMAX)
                        ++count;
                    //printf("row=%lu, count=%lu\n", i, count);
                    //printf("\tFQHE:\t%d %d %d %d \ttype:\t%d\n", j1 ,j2, j3, j4, layerQ(j1, j2, (int)N_phi));
                }
            }
        }
    }
    //printf("row=%lu, count=%lu\n", i, count);
    free(bra_i);
    free(bra_j);
    free(diff);
    return count+1;
}


size_t tight_bond_Z2_count(const Data_basis *basis, Data_basis *bra_i, Data_basis *bra_j, const hash & hashT, size_t N_o, size_t N_e)
{
    size_t j1, j2, count = 0;
    size_t N_orbit = 2*N_o;
    memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
    memset(bra_j, 0, (size_t)N_e * sizeof(Data_basis));
    for(size_t i=0;i<N_e;++i)
        bra_i[basis[i]] = 1;
    int n_find;

    //printf("dict_num=%lu, binary=%lu\n", dict_num(basis, N_e), dict_num_binary(bra_i, N_e, N_orbit));
    for(size_t i=0;i<N_e;++i)
    {
        j2 = basis[i];
        j1 = tight_up(j2, N_o);
        if(bra_i[j1]==0)
        {
            bra_i[j2] = 0;
            bra_i[j1] = 1;

            binary_to_basis(bra_i, bra_j, N_orbit);
            auto idx = find_Z2_index(bra_j, hashT, N_e, N_o,n_find);

            /*for(size_t i=0;i<N_e;++i)
                printf("%u ", basis[i]);
            printf("\t%lu->%lu\t", j2, j1);
            for(size_t i=0;i<N_e;++i)
                printf("%u ", bra_j[i]);*/
            if(idx!=ULIMAX)
            {
                ++count;
                //printf("dict_num=%lu, binary=%lu, Z2_index=%lu\n", dict_num(basis, N_e), dict_num_binary(bra_i, N_e, N_orbit), idx);
            }
            //else
                //printf("dict_num=%lu, binary=%lu, NOT FOUND\n", dict_num(basis, N_e), dict_num_binary(bra_i, N_e, N_orbit));

            bra_i[j2] = 1;
            bra_i[j1] = 0;
        }
    }

    return count;
}
/*
size_t tight_up(size_t j1, size_t N_phi)
{
    if(j1>=N_phi)
        return j1 - N_phi;
    else
        return j1 + N_phi;
}*/



extern "C" size_t twoBodyZ2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double t, int num_thread)
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
        auto count = twoBodyZ2Generate_thread(indptr, indices, data, A1_list, A12_list, \
                                    A.basis, A.hashT, A.Z2_phase, A.Z_2, t, A.dim, A.N_o, A.N_e, i);
        if(count!= indptr[i+1] - indptr[i])
        {
            printf("In function twoBodyZ2Generate: count=%lu, but indptr[%lu] - indptr[%lu] = %lu\n", count, i+1,i,indptr[i+1] - indptr[i]);
            exit(10);
        }
    }

    return indptr[A.dim];
}

size_t twoBodyZ2Generate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, \
                    Data_basis *k1_basis, const hash & hashT, char* Z2_phase, int Z_2, double t, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row)
{
    size_t bra_num, ind, count = 0;
    Data_basis *basis_t = k1_basis + row * N_e;
    size_t j1, j2, j3, j4, j13, j14, j23, j24, J12, J34, J1, J2, J3, J4;
    size_t N_orbit = 2 * N_o;
    Data_dtype dig = 0.0;
    MKL_INT *indices_thread = indices + indptr[row];
    Data_dtype* data_thread = data + indptr[row];
    int phase, n_find;

    Data_basis *bra_i = (Data_basis *)malloc(N_orbit * sizeof(Data_basis));
    Data_basis *bra_j = (Data_basis *)malloc(N_e * sizeof(Data_basis));
    Data_basis *diff = (Data_basis *)malloc((N_orbit-N_e) * sizeof(Data_basis));
    // CDW term
    //cdw_energy(k1_basis+i*N_e, N_phi, N_e); //dig
    //printf("\t%u %lu\n", index[dict_num(basis_t, N_e)], i);
    
    // tight bond term
    count += tight_bond_Z2(basis_t, bra_i, bra_j, hashT, Z2_phase, Z_2, indices_thread, data_thread, \
                        t, N_o, N_e, row);// update count
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
                    intra_hopping(basis_t, bra_i, bra_j, 2*N_o, N_e, j1, j2, j3, j4, &phase);
                    ind = find_Z2_index(bra_j, hashT, N_e, N_o, n_find);
                    if(ind==ULIMAX)
                        continue;

                    J1 = j1%N_o;
                    J2 = j2%N_o;
                    J3 = j3%N_o;
                    J4 = j4%N_o;
                    if( n_find==2 && Z_2*Z2_phase[ind]<0 )
                        phase *= -1;
                    if(layerQ(j1, j2, N_o)==0)
                        dig = 2.0 * A12_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4] * phase;
                    else
                        dig = 2.0 * (A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4]-A1_list[J2*N_o*N_o*N_o + J1*N_o*N_o + J3*N_o + J4]) * phase;
                    dig *= std::sqrt( std::abs((double)Z2_phase[row]/Z2_phase[ind]) );
                    //printf("off:\t\t\b (%d, %d, %d, %d),(%d, %d, %d, %d), ind=%lu row=%lu layerQ=%d: ", \
                                j1, j2, j3, j4,J1,J2,J3,J4, ind, row, layerQ(j1, j2, N_o));
                    //printf("%lf=%lf*%lf, phase=%d\n", dig, A12_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4], std::sqrt( std::abs((double)Z2_phase[row]/Z2_phase[ind]) ), phase);
                    //printf("%lf, %lf, %lf, %lu\n", A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4], A1_list[J2*N_o*N_o*N_o + J1*N_o*N_o + J3*N_o + J4], A1_list[17], u);
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
            if(layerQ(j1, j2, N_o)==0)
                dig += 2.0 * A12_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4];
            else
                dig += 2.0 * (A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4]-A1_list[J2*N_o*N_o*N_o + J1*N_o*N_o + J3*N_o + J4]);
            //printf("diag: %d, %d, %d, %d, %lu\n", j1, j2, j3, j4, row);
            //printf("diag: %lf+I%lf\n", dig, 0);
            //printf("diag: %d %d %d\n",layerQ(j3, j4, (int)N_o),((j2<2*N_o) && (j2>=N_o)), j2>=N_o);
                    //printf("diag:\t\t\b (%d, %d, %d, %d), row=%lu layerQ=%d: ", \
                                j1, j2, j3, j4, row, layerQ(j1, j2, N_o));
                    //printf("%lf+I%lf, %lf\n", dig, 0.0, A1_list[J1*N_o*N_o*N_o + J2*N_o*N_o + J3*N_o + J4]);
        }
    }
    //printf("%lu %lu %u %lf+I%lf\n", row, count, indptr[row+1], real(dig), imag(dig));
    //getchar();
    indices_thread[count] = row;
    data_thread[count] = dig;// + cdw_energy(basis_t, cdw, N_o, N_e);
    //printf("%lu %lf %lf+I%lf\n", row, cdw_energy(basis_t, cdw, N_phi, N_e), creal(dig), cimag(dig));
    ++count;


    free(bra_i);
    free(bra_j);
    free(diff);
    return count;
}

size_t tight_bond_Z2(const Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
                const hash & hashT, const char* Z2_phase, int Z_2, MKL_INT *indices_thread, Data_dtype *data_thread, \
                double t, size_t N_o, size_t N_e, size_t row)
{
    size_t j1, j2, count = 0;
    size_t N_orbit = 2*N_o;
    memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));
    memset(bra_j, 0, (size_t)N_e * sizeof(Data_basis));
    for(size_t i=0;i<N_e;++i)
        bra_i[basis_t[i]] = 1;
    int phase, b, e, n_find;

    //size_t ket_num = dict_num(basis_t, N_e);
    //printf("dict_num=%lu, binary=%lu\n", dict_num(basis, N_e), dict_num_binary(bra_i, N_e, N_orbit));
    for(size_t i=0;i<N_e;++i)
    {
        j2 = basis_t[i];
        j1 = tight_up(j2, N_o);
        if(bra_i[j1]==0)
        {
            bra_i[j2] = 0;
            bra_i[j1] = 1;

            binary_to_basis(bra_i, bra_j, N_orbit);
            auto idx = find_Z2_index(bra_j, hashT, N_e, N_o, n_find);

            //printf("hopping:\t");
            //for(size_t i=0;i<N_e;++i)
            //    printf("%u ", basis_t[i]);
            //printf("\t%lu->%lu\t", j2, j1);
            //for(size_t i=0;i<N_e;++i)
            //    printf("%u ", bra_j[i]);
            if(idx!=ULIMAX)
            {
                //size_t bra_num = dict_num(basis_t, N_e);
                b = (j1>j2)?j2:j1;
                e = (j1<j2)?j2:j1;
                phase = 0;
                for(size_t cc=b+1;cc<e;++cc)
                    phase += bra_i[cc];
                if( n_find==2 && Z_2*Z2_phase[idx]<0 )
                    phase += 1;
                phase = (phase&1)?-1:1;

                indices_thread[count] = idx;
                data_thread[count] = t * phase * std::sqrt( std::abs((double)Z2_phase[row]/Z2_phase[idx]) );
                //printf("(dict_num=%lu,row=%lu)->(dict_num=%lu,ind=%lu,n_find=%d),count=%lu, data=%lf\n", dict_num(basis_t, N_e),row, dict_num_binary(bra_i, N_e, N_orbit), idx,n_find,count,data_thread[count]);
                ++count;
            }
            //else
                //printf("dict_num=%lu, binary=%lu, NOT FOUND\n", dict_num(basis, N_e), dict_num_binary(bra_i, N_e, N_orbit));

            bra_i[j2] = 1;
            bra_i[j1] = 0;
        }
    }

    return count;
}