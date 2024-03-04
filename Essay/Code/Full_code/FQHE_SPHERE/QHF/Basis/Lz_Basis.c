

int fermionLzBasis(size_t N_phi, size_t N_e, int L_z, char *temp,\
					char *basis, const size_t N_p, size_t *o)
{
	if(N_e==0)
	{
		if(L_z!=0)
			return 0;
		char *basis_t = basis + (size_t)o[0] * N_p;
		for(size_t j=0;j<N_p;++j)
		{
			basis_t[j] = temp[j];
			//printf("%d ",temp[j]);
		}
		//printf("\t%d\n", L_z);
		++o[0];
		return 0;
	}

	if(L_z<0)
		return 0;
	size_t j;
	for(size_t i=N_e-1;i<N_phi;++i)
	{
		j = (i>=N_p)?(i-N_p):i;
		temp[(size_t)N_e - 1] = i;
		fermionLzBasis(i, N_e-1, L_z-j, temp, basis, N_p, o);
	}
	return 0;
}

int fermionLzBasisCount(size_t N_phi, size_t N_e, int L_z, char *temp,\
					const size_t N_p, size_t *o)
{
	if(N_e==0)
	{
		if(L_z!=0)
			return 0;
		++o[0];
		return 0;
	}

	if(L_z<0)
		return 0;
	size_t j;
	for(size_t i=N_e-1;i<N_phi;++i)
	{
		j = (i>=N_p)?(i-N_p):i;
		fermionLzBasisCount(i, N_e-1, L_z-j, temp, N_p, o);
	}
	return 0;
}


size_t dict_num(char *basis, size_t N_e)
{
	size_t num = 0;
	for(size_t i=0;i<N_e;++i) // particle number
		num += binomial(basis[i], i+1);
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

size_t *basis_to_num(char *basis, size_t dim, size_t N_e, MKL_INT *index)
{
	char *basis_t = basis;
	size_t *base_num = NULL;

	base_num = (size_t *)malloc(dim * sizeof(size_t));
	if(base_num==NULL)
		return NULL;

	for(size_t i=0;i<dim;++i)
	{
		base_num[i] = dict_num(basis_t, N_e);
		index[(size_t)base_num[i]] = i;
		basis_t += (size_t)N_e;
	}
	return base_num;
}