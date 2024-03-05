

//void LzZ2_to_Lz(const Data_dtype* LzZ2_ket, Data_dtype* Lz_ket, size_t )

// you should create B struct first, 
// init_FQHE_B(size_t N_o, size_t N_e);
// create_Lz_B(size_t L_z);
extern "C" void LzZ2PH_to_Lz(const Data_dtype* LzZ2PH_ket, Data_dtype* Lz_ket, int num_thread)
{
    # pragma omp parallel for num_threads(num_thread)
    for(size_t i=0;i<get_B_Lz_dim();++i)
        LzZ2PH_to_Lz_thread(LzZ2PH_ket, Lz_ket, i);
    return;
}


void LzZ2PH_to_Lz_thread(const Data_dtype* LzZ2PH_ket, Data_dtype* Lz_ket, size_t row)
{
    int length = 1;
    char phase = 1;
    size_t ind = find_Z2PH_index(B.basis+(size_t)B.N_e*row, A.hashT, A.PHhashT, A.Z2_phase, A.PH_phase, 
                                A.Z2PH_length, A.PHindex, A.Z_2, A.PH, A.N_e, A.N_o, length, phase);
    if(ind!=ULIMAX)
    {
        Lz_ket[row] = LzZ2PH_ket[ind] * phase/std::sqrt( std::abs(length) );
    }
    return;
}