

// c_up^dagger -> c_down, c_down^dagger -> -c_up
void PH_transformation(const Data_basis *basis, Data_basis *PHbasis, size_t N_o, size_t N_e, char *phase)
{
    Data_basis *temp =  (Data_basis *)malloc((size_t)N_e * sizeof(Data_basis));
    Z2_transformation(basis, temp, N_o, N_e, phase);//*phase = ((h*(N_e-h))&1==1)?-1:1;

    for(size_t i=0;i<N_e;++i)
        *phase *= (basis[i]<N_o)?1:-1;//c_down^dagger -> -c_up


    int new_phase = 0;
    size_t di=0, dj=0, N_orbit=2*N_o;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(temp[di]==i)
        {
            ++di;
            continue;
        }
        PHbasis[dj] = i;
        new_phase += i;
        ++dj;
    }
    *phase *= (new_phase&1==1)?-1:1;

    free(temp);
}

void Z2_transformation(const Data_basis *basis, Data_basis *Z2basis, size_t N_o, size_t N_e, char *phase)
{
    size_t p_j, j, h, N_orbit=2*N_o;
    size_t num = 0;
    for(h=0;h<N_e;++h)
    {
        if(basis[h]+N_o>=N_orbit)
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
        p_j = basis[j] + N_o;
        if(p_j>=N_orbit)
            p_j -= N_orbit;
        Z2basis[i] = p_j;
        //num += binomial(p_j, i+1);
    }
    return;
}


// the phase may not be exact, another function with Z_2 arguement is exact
size_t find_Z2_index(const Data_basis *basis, const hash & hashT, const char* Z2_phase, size_t N_e, size_t N_o, int& length, char& phase)
{
    int n_find = 1;
    size_t num = dict_num(basis, N_e);
    auto it = hashT.find( num );
    if(it!=hashT.end())
    {
        phase = 1;
        length = std::abs(Z2_phase[it->second]);
        return it->second;
    }
    n_find = 2;
    num = base_translation(basis, N_o, 2*N_o, N_e, &phase);
    it = hashT.find( num );
    if(it!=hashT.end())
    {
        length = std::abs(Z2_phase[it->second]);
        return it->second;
    }
    else
        return ULIMAX;
}


// the phase is exact
size_t find_Z2_index(const Data_basis *basis, const hash & hashT, const char* Z2_phase, size_t N_e, size_t N_o, int Z_2, int& length, char& phase)
{
    int n_find = 1;
    size_t num = dict_num(basis, N_e);
    auto it = hashT.find( num );
    if(it!=hashT.end())
    {
        phase = 1;
        length = std::abs(Z2_phase[it->second]);
        return it->second;
    }
    n_find = 2;
    num = base_translation(basis, N_o, 2*N_o, N_e, &phase);
    it = hashT.find( num );
    if(it!=hashT.end())
    {
        //if( Z2_phase[it->second]*Z_2<0 )
        phase *= Z_2;
        length = std::abs(Z2_phase[it->second]);
        return it->second;
    }
    else
        return ULIMAX;
}

// the phase may not be exact, another function with Z_2 and PH arguements is exact
size_t find_Z2PH_index(const Data_basis *basis, const hash & hashT, const hash & PHhashT, const char* Z2_phase,\
                    const char* PH_phase, const char* Z2PH_length, const size_t* PHindex, size_t N_e, size_t N_o, int& length, char& phase)
{
    phase = 1;
    char tphase = 1;
    int n_find = 1;
    length = 1;
    size_t num = find_Z2_index(basis, hashT, Z2_phase, N_e, N_o, length, tphase);
    auto it = PHhashT.find(num);
    if( num!=ULIMAX && it!=PHhashT.end() )
    {
        phase = tphase;
        length *= PH_phase[it->second];
        return it->second;
    }

    Data_basis *PHbasis =  (Data_basis *)malloc((size_t)N_e * sizeof(Data_basis));
    PH_transformation(basis, PHbasis, N_o, N_e, &phase);
    num = find_Z2_index(PHbasis, hashT, Z2_phase, N_e, N_o, length, tphase);
    free(PHbasis);
    it = PHhashT.find(num);
    if( num!=ULIMAX && it!=PHhashT.end() )
    {
        phase *= tphase;
        length *= PH_phase[it->second];
        return it->second;
    }

    return ULIMAX;
}

// the phase is exact
size_t find_Z2PH_index(const Data_basis *basis, const hash & hashT, const hash & PHhashT, const char* Z2_phase,\
                    const char* PH_phase, const char* Z2PH_length, const size_t* PHindex, int Z_2, int PH, \
                        size_t N_e, size_t N_o, int& length, char& phase)
{
    phase = 1;
    char tphase = 1;
    int n_find = 1;
    length = 1;
    size_t num = find_Z2_index(basis, hashT, Z2_phase, N_e, N_o, Z_2, length, tphase);
    auto it = PHhashT.find(num);
    if( num!=ULIMAX && it!=PHhashT.end() )
    {
        phase = tphase;
        length *= PH_phase[it->second];
        return it->second;
    }

    Data_basis *PHbasis =  (Data_basis *)malloc((size_t)N_e * sizeof(Data_basis));
    PH_transformation(basis, PHbasis, N_o, N_e, &phase);
    num = find_Z2_index(PHbasis, hashT, Z2_phase, N_e, N_o, Z_2, length, tphase);
    free(PHbasis);
    it = PHhashT.find(num);
    if( num!=ULIMAX && it!=PHhashT.end() )
    {
        phase = phase * tphase;
        //if( PH_phase[it->second]*PH<0 )
        phase *= PH;
        length *= PH_phase[it->second];
        return it->second;
    }

    return ULIMAX;
}

