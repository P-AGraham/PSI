
extern "C" size_t binomial(size_t m, size_t n)//C_m^n
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

extern "C" int fermion_basis(size_t N_phi, size_t N_e, unsigned char *temp,\
					unsigned char *basis, const size_t N_p, size_t *o)
{
	if(N_e==0)
	{
		unsigned char *basis_t = basis + (size_t)o[0] * N_p;
		for(size_t j=0;j<N_p;++j)
		{
			basis_t[j] = temp[j];
			//printf("%d ",temp[j]);
		}
		//printf("\t%lu\n",1);//dict_num(temp,N_p));
		++o[0];
		return 0;
	}

	for(size_t i=N_e-1;i<N_phi;++i)
	{
		temp[(size_t)N_e - 1] = i;
		fermion_basis(i, N_e-1, temp, basis, N_p, o);
	}
	return 0;
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

size_t dict_num(const unsigned char* basis, size_t N_e)
{
	size_t num = 0;
	for(size_t i=0;i<N_e;++i) // particle number
		num += binomial(basis[i], i+1);
	return num;
}

size_t dict_num_minus_phi(const unsigned char* basis, size_t N_e, size_t N_phi)
{
	size_t num = 0;
	for(size_t i=0;i<N_e;++i) // particle number
		num += binomial(basis[i]-N_phi, i+1);
	return num;
}

size_t base_translation(const unsigned char* basis, size_t offset,\
						size_t N_phi, size_t N_e, char *phase)
{
	size_t p_j, j, h;
	size_t num = 0;
	for(h=0;h<N_e;++h)
	{
		if(basis[h]+offset>=N_phi)
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
		if(p_j>=N_phi)
			p_j -= N_phi;
		num += binomial(p_j, i+1);
	}
	return num;
}