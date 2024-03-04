
size_t hopping(const unsigned char *basis, unsigned char *bra_i, unsigned char *bra_j, \
					size_t N_phi, size_t N_e, int j1, int j2, int j3, int j4, int *phase)
{
	// j4 > j3, j1 > j2

	// initialization
	memset(bra_i, 0, (size_t)N_phi * sizeof(unsigned char));
	//for(size_t i=0;i<N_phi;++i)
	//	bra_i[i] = 0;

	// translation occupy -> binary
	for(size_t i=0;i<N_e;++i)
		++bra_i[basis[i]];
	*phase = 0;

	// annihilation
	for(size_t i=j3+1;i<j4;++i)
		*phase += bra_i[i];
	--bra_i[j4];
	--bra_i[j3];

	// creation
	for(size_t i=j2+1;i<j1;++i)
		*phase += bra_i[i];
	++bra_i[j2];
	++bra_i[j1];

	*phase = (*phase&1)?-1:1;
	
	// translation binary -> occupy
	size_t k = 0;
	for(size_t i=0;i<N_phi;++i)
	{
		if(bra_i[i])
		{
			bra_j[k] = i;
			++k;
		}
	}

	return dict_num(bra_j, N_e);
}

size_t hoppingOneBody(const unsigned char *basis, unsigned char *bra_i, unsigned char *bra_j, \
					size_t N_phi, size_t N_e, int j1, int j2, int *phase)
{
	// j4 > j3, j1 > j2

	// initialization
	memset(bra_i, 0, (size_t)N_phi * sizeof(unsigned char));
	//for(size_t i=0;i<N_phi;++i)
	//	bra_i[i] = 0;

	// translation occupy -> binary
	for(size_t i=0;i<N_e;++i)
		++bra_i[basis[i]];
	*phase = 0;

	// annihilation
	--bra_i[j2];
	// creation
	++bra_i[j1]; 

	size_t left, right;
	if(j1>j2)
	{
		left = j2;
		right = j1;
	}
	else
	{
		left = j1;
		right = j2;
	}
	for(size_t i=left+1;i<right;++i)
		*phase += bra_i[i];

	*phase = (*phase&1)?-1:1;
	
	// translation binary -> occupy
	size_t k = 0;
	for(size_t i=0;i<N_phi;++i)
	{
		if(bra_i[i])
		{
			bra_j[k] = i;
			++k;
		}
	}

	return dict_num(bra_j, N_e);
}
