
void DLFermion::k2_to_k1_transformation(complex<double>* ket_k2, complex<double>* ket_k1)
{
	# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<A.get_k1_dim();++i)
		H_space.k2_to_k1_thread(i, ket_k2, ket_k1);
	return;
}

void basis_symm::k2_to_k1_thread(size_t i, complex<double>* ket_k2, complex<double>* ket_k1)
{
	unordered_map<size_t, size_t>::iterator hash_ed = k2_hash.end();
	size_t k1_ind = k1_basis[i].first;
	unordered_map<size_t, size_t>::iterator hash_it = k2_hash.find(k1_ind);

	if(hash_it!=hash_ed) // filter
	{
		size_t k2_ind = hash_it->second;
		ket_k1[i] = ket_k2[k2_ind] * exp(2.0*M_PI*ii*k2*k1_basis[i].sn/(pqgcd)) / sqrt((double)k1_basis[k1_ind].period);
		if(k1_basis[i].period==0)
			ket_k1[i] *= k1_basis[i].phase;
	}
}

void basis_symm::occupation_number(complex<double>* ket_k1, double *out, int num_thread)
{
	double temp;
	for(size_t i=0;i<k1_dim;++i)
	{
		this->set(i);
		temp = norm(ket_k1[i]);
		for(size_t j=0;j<N_u;++j)
			out[ket->up()[j]] += temp;
		for(size_t j=0;j<N_d;++j)
			out[N_phi+ket->down()[j]] += temp;
	}
	this->set(0);
}

size_t basis_symm::to_3layers_trans(complex<double>* ket_k1, complex<double> *out, int num_thread)
{
	basis_array tlayer_basis = basis_array();
	tlayer_basis.create(N_phi*3, N_u+N_d, num_thread, N_phi);
	size_t key, ind, cc=0;
	for(size_t i=0;i<tlayer_basis.get_dim();++i)
	{
		if(tlayer_basis.get_k(i) == k1)
		{
			if(tlsyers_Q(tlayer_basis[i], N_phi, N_u, N_d))
			{
				key = dict_num(tlayer_basis[i], N_u)*stride + dict_num_minus_phi(tlayer_basis[i]+N_u, N_d, N_phi);
				ind = find(key);
				out[cc] += ket_k1[ind];
			}
			++cc;
		}
	}
	return cc;
}

bool tlsyers_Q(const unsigned char *basis, size_t N_phi, size_t N_u, size_t N_d)
{
	size_t c_up, c_down;
	c_up = c_down = 0;
	for(size_t i=0;i<N_u+N_d;++i)
	{
		//cout << (int)basis[i] << " ";
		if(basis[i]<N_phi)
			++c_up;
		else if(basis[i]<N_phi*2)
			++c_down;
	}
	//cout << "N_up = " << c_up << " N_down = " << c_down << "\t";
	//cout << ((c_up==N_u) && (c_down==N_d));
	return (c_up==N_u) && (c_down==N_d);
}