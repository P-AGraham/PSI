

size_t base_translation(const Data_basis *basis, size_t offset, size_t N_orbit, size_t N_e, char *phase)
{
    size_t p_j, j, h;
    size_t num = 0;
    for(h=0;h<N_e;++h)
    {
        if(basis[h]+offset>=N_orbit)
            break;
    }
    if(h>=N_e)
        h -= N_e;
    *phase = ((h*(N_e-h))&1==1)?-1:1;
    for(size_t i=0;i<N_e;++i)
    {
        j = i + h;
        if(j>=N_e)
            j -= N_e;
        p_j = basis[j] + offset;
        if(p_j>=N_orbit)
            p_j -= N_orbit;
        num += binomial(p_j, i+1);
    }
    return num;
}

int fermionLzZ2Basis(Data_basis N_tot_o, Data_basis N_e, int Lz, const int Z2, Data_basis *temp,\
                    Data_basis *basis, char* Z2_phase, size_t *Z2_num, const Data_basis N_p, const Data_basis N_o, bool projectionQ, size_t *o)//N_tot_o = 2*N_o
{
    if(N_e==0)
    {
        if(Lz!=0)
            return 0;
        char phase;
        size_t bra_num = base_translation(temp, N_o, 2*N_o, N_p, &phase);
        size_t ket_num = dict_num(temp,N_p);
        if(bra_num<ket_num) return 0;
        int state = phase*( (bra_num==ket_num)?1:2 );

        if(Z2==1)// P(+)|phase=+1> and P(-)|phase=-1, T=2>
        {
            if(state==-1) 
                return 0;
        }
        else if(Z2==-1)// P(-)|phase=+1> and P(+)|phase=-1, T=2> and P(-)|phase=+1, T=2>
        {
            if(state==1) 
                return 0;
        }
        else
            return 0;

        if( projectionQ && !projectionTest(temp, N_p, N_o) )
            return 0;

        Z2_phase[(*o)] = state;
        Data_basis *basis_t = basis + (size_t)o[0] * N_p;
        for(size_t j=0;j<N_p;++j)
        {
            basis_t[j] = temp[j];
            //printf("%d ",temp[j]);
        }//printf(": %+d, %d, %d\n", state, projectionQ, projectionQ && !projectionTest(temp, N_p, N_o));
        #ifdef __PRINT_Z2INDEX
        Z2_num[*o] = bra_num;
        #endif
        ++o[0];
        return 0;
    }

    if(Lz<0)
        return 0;
    for(Data_basis i=N_e-1;i<N_tot_o;++i)
    {
        temp[(size_t)N_e - 1] = i;
        fermionLzZ2Basis(i, N_e-1, Lz-(i%N_o), Z2, temp, basis, Z2_phase, Z2_num, N_p, N_o, projectionQ, o);
    }
    return 0;
}

int fermionLzZ2BasisCount(Data_basis N_tot_o, Data_basis N_e, int Lz, const int Z2, Data_basis *temp,\
                      const Data_basis N_p, const Data_basis N_o, bool projectionQ, size_t *o)//N_tot_o = 2*N_o
{
    if(N_e==0)
    {
        if(Lz!=0)
            return 0;
        char phase;
        size_t bra_num = base_translation(temp, N_o, 2*N_o, N_p, &phase);
        size_t ket_num = dict_num(temp,N_p);
        if(bra_num>ket_num) return 0;
        int state = phase*( (bra_num==ket_num)?1:2 );

        if(Z2==1)// P(+)|phase=+1> and P(-)|phase=-1, T=2>
        {
            if(state==-1) 
                return 0;
        }
        else if(Z2==-1)// P(-)|phase=+1> and P(+)|phase=-1, T=2> and P(-)|phase=+1, T=2>
        {
            if(state==1) 
                return 0;
        }
        else
            return 0;

        if( projectionQ && !projectionTest(temp, N_p, N_o) )
            return 0;

        //for(size_t j=0;j<N_p;++j)
        //{
        //    printf("%d ",temp[j]);
        //}printf("\n");
        ++o[0];
        return 0;
    }

    if(Lz<0)
        return 0;
    for(Data_basis i=N_e-1;i<N_tot_o;++i)
    {
        temp[(size_t)N_e - 1] = i;
        fermionLzZ2BasisCount(i, N_e-1, Lz-(i%N_o), Z2, temp, N_p, N_o, projectionQ, o);
    }
    return 0;
}


size_t find_Z2_index(const Data_basis *basis, const hash & hashT, size_t N_e, size_t N_o, int& n_find)
{
    n_find = 1;
    size_t num = dict_num(basis, N_e);
    auto it = hashT.find( num );
    if(it!=hashT.end())
        return it->second;
    n_find = 2;
    char phase;
    num = base_translation(basis, N_o, 2*N_o, N_e, &phase);
    it = hashT.find( num );
    if(it!=hashT.end())
        return it->second;
    else
        return ULIMAX;
}