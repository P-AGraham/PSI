
//--------------------------------------------------------------------------------------------------------------------
//-------------------------<class basis_array>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

size_t basis_array::create(size_t N_orbit, size_t N_ele, int num_thread, size_t N_mod)
{
	if(basis!=NULL)
		delete basis;
	if(moment!=NULL)
		delete moment;
	N_e = N_ele;
	N_phi = N_orbit;
	dim = binomial(N_phi, N_e);
	basis = new unsigned char [dim*N_e];
	unsigned char *temp = new unsigned char [N_e];
	size_t o = 0;
	fermion_basis(N_phi, N_e, temp, basis, N_e, &o);
	delete [] temp;

	if(N_mod==0)
		N_mod = N_phi;

	moment = new unsigned char [dim];
	# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<dim;++i)
		moment[i] = cal_moment_thread(basis+i*N_e, N_e, N_phi, N_mod);
	return dim;
}

std::ostream & operator<<(std::ostream & os, const basis_array & b)
{
	for(size_t i=0;i<b.dim;++i)
	{
		for(size_t j=0;j<b.N_e;++j)
			os << (int)b.basis[i*b.N_e+j] << " ";
		os << "\tk = " << (int)b.moment[i] << "\n";
	}
	os << "dim = " << b.dim;
	return os;
}

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<function>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------


inline unsigned char cal_moment_thread(unsigned char *basis, unsigned char N_e, unsigned char N_phi, unsigned char N_m)
{
	int m = 0;
	for(size_t i=0;i<N_e;++i)
		m += basis[i];
	return (unsigned char)(m%N_m);
}