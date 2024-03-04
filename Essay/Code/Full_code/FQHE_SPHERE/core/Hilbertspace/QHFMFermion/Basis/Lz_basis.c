//  3Body_Fermion_k1_basis.h

int fermionLzBasis(Data_basis N_tot_o, Data_basis N_e, int Lz, Data_basis *temp,\
					Data_basis *basis, const Data_basis N_p, const Data_basis N_o, size_t *o)//N_tot_o = 2*N_o
{
	if(N_e==0)
	{
		if(Lz!=0)
			return 0;
		Data_basis *basis_t = basis + (size_t)o[0] * N_p;
		for(size_t j=0;j<N_p;++j)
		{
			basis_t[j] = temp[j];
			//printf("%d ",temp[j]);
		}
		//printf("\t%lu\n",dict_num(temp,N_p));
		++o[0];
		//getchar();
		return 0;
	}

	if(Lz<0)
		return 0;
	for(Data_basis i=N_e-1;i<N_tot_o;++i)
	{
		temp[(size_t)N_e - 1] = i;
		fermionLzBasis(i, N_e-1, Lz-(i%N_o), temp, basis, N_p, N_o, o);
	}
	return 0;
}


int fermionLzBasisCount(Data_basis N_tot_o, Data_basis N_e, int Lz, const Data_basis N_p, const Data_basis N_o, size_t *o)
{
	if(N_e==0)
	{
		if(Lz!=0)
			return 0;
		++o[0];
		return 0;
	}

	if(Lz<0)
		return 0;
	for(Data_basis i=N_e-1;i<N_tot_o;++i)
		fermionLzBasisCount(i, N_e-1, Lz-(i%N_o), N_p, N_o, o);
	return 0;
}

size_t binomial(size_t m, size_t n)//C_m^n
{
	size_t b;
	if(m<n||n<0)
		return 0;
	b = Binomial[m][n];
	if(b)
		return b;
	if(m<n)
		return 0;
	if(n==0)
		return 1;
	if(n==1)
		return m;
	if(n>m/2)
	{
		b = binomial(m, m-n);
		Binomial[m][n] = b;
		return b;
	}
	else
	{
		b = binomial(m-1, n-1)+binomial(m-1, n);
		Binomial[m][n] = b;
		return b;
	}
}

size_t dict_num(const Data_basis *basis, size_t N_e)
{
	size_t num = 0;
	for(size_t i=0;i<N_e;++i) // particle number
		num += binomial(basis[i], i+1);
	return num;
}

size_t dict_num_binary(const Data_basis *basis, size_t N_e, size_t N_o_tot)
{
	size_t num = 0, i=0;
	for(size_t j=0;j<N_o_tot;++j)
	{
		if(basis[j])
		{
			++i;
			num += binomial(j, i);
			//printf("binomial(%lu, %lu) + ", j, i);
		}
	}
		//printf("\n");
	return num;
}

size_t gcd(size_t x, size_t y)
{
	size_t a, b;
	a = x;
	b = y;
	while(a != b)
	{
		if(a > b)
			a -= b;
		else
			b -= a;
	}
	return a;
}


void basis_to_num(Data_basis *basis, size_t dim, size_t N_e, std::unordered_map<size_t, size_t> & hashT)
{
	Data_basis *basis_t = basis;
	size_t base_num;

	for(size_t i=0;i<dim;++i)
	{
		base_num = dict_num(basis_t, N_e);
		hashT.insert({base_num, i});
		basis_t += (size_t)N_e;
	}
	return;
}

void binary_to_basis(Data_basis *ket_binary, Data_basis *ket_basis, size_t N_orbit)
{
    // translate binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(ket_binary[i])
        {
            ket_basis[k] = i;
            ++k;
        }
    }
}

void basis_to_binary(const Data_basis *ket_basis, Data_basis *ket_binary, size_t N_orbit, size_t N_e)
{
    // translate occupy -> binary
    memset(ket_binary, 0, (size_t)N_orbit * sizeof(Data_basis));
    for(auto i=0;i<N_e;++i)
    	ket_binary[ket_basis[i]] = 1;
}

void print_basis(const Data_basis *basis, size_t N_e, size_t N_o)
{
	for(size_t j=0;j<N_e;++j)
	{
		if(basis[j]<N_o)
			printf("%d ",basis[j]);
		else
			printf("\033[1;31m%d \033[0m",basis[j]);
	}
}

void println_basis(const Data_basis *basis, size_t N_e, size_t N_o)
{
	print_basis(basis, N_e, N_o);
	printf("\n");
}

void print_basis_mod(const Data_basis *basis, size_t N_e, size_t N_o)
{
	for(size_t j=0;j<N_e;++j)
	{
		if(basis[j]<N_o)
			printf("%d ",basis[j]);
		else
			printf("\033[1;31m%d \033[0m",basis[j]%N_o);
	}
}

void println_basis_mod(const Data_basis *basis, size_t N_e, size_t N_o)
{
	print_basis_mod(basis, N_e, N_o);
	printf("\n");
}

bool projectionTest(const Data_basis *basis, size_t N_e, size_t N_o)
{
	size_t N_orbit = 2*N_o;
	Data_basis *bin_basis = (Data_basis *)malloc((size_t)N_orbit * sizeof(Data_basis));
	basis_to_binary(basis, bin_basis, N_orbit, N_e);
	for(auto i=0;i<N_o;++i)
	{
		if(bin_basis[i]==bin_basis[i+N_o])
		{
			free(bin_basis);
			return false;
		}
	}

	free(bin_basis);
	return true;
}