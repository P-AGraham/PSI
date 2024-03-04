

complex<double> ii = {0, 1};

class DLFermion
{
private:
	size_t N_phi_;
	size_t N_e_;
	size_t Lz_ = 0xFFFFFFFFFFFFFFFF;
	int num_thread_;

	basis_array basis_ = basis_array();

	//basis_symm H_space;
public:
	DLFermion(size_t N_orbit, size_t N_e, int num_threads);
	~DLFermion();
	size_t get_N_phi(){return N_phi_;}
	size_t get_N_e(){return N_e_;}
	//size_t get_Lz_dim(){return H_space.get_Lz_dim();}

	void reset(size_t N_orbit, size_t N_e, int num_threads);
	basis_array & get_basis_1(){return basis_;}
	//basis_symm & get_basis_symm(){return H_space;}

	//void create_Lz_basis(size_t Lz);
	//size_t count_2symm(unsigned int *indptr){return H_space.element_count(indptr, num_thread);}
	//size_t generate_2symm(complex<double>* A1_list, complex<double>* A2_list, complex<double>* A12_list, \
						unsigned int* indptr, unsigned int* indices, complex<double>* data)
	//{return H_space.element_generate(A1_list, A2_list, A12_list, indptr, indices, data, num_thread);}

	//void copy_basis(unsigned char *out, size_t i){H_space.copy_basis(out, i);}
	//void occupation_number(complex<double>* ket, double *out)
	//{H_space.occupation_number(ket, out, num_thread);}
	//size_t to_3layers_trans(complex<double>* ket_k1, complex<double> *out){return H_space.to_3layers_trans(ket_k1, out, num_thread);}

	friend std::ostream & operator<<(std::ostream & os, const DLFermion & b);
};

DLFermion::DLFermion(size_t N_orbit, size_t N_ele, int num_threads)
{
	N_phi_ = N_orbit;
	N_e_ = N_e;
	num_thread_ = num_threads;
	Lz_ = 0xFFFFFFFFFFFFFFFF;
	basis_ = basis_array();
	std::cout << "create DLFermion at " << this << std::endl;
}

void DLFermion::reset(size_t N_orbit, size_t N_up, size_t N_down, int num_threads)
{
	N_phi_ = N_orbit;
	N_e_ = N_e;
	num_thread_ = num_threads;
	Lz_ = 0xFFFFFFFFFFFFFFFF;
	basis_.reset();
	//H_space.reset();
}

DLFermion::~DLFermion()
{
	std::cout << "destroy DLFermion at " << this << std::endl;
}

std::ostream & operator<<(std::ostream & os, const DLFermion & b)
{
	os << "N_phi = " << b.N_phi_ << " N_e = " << b.N_e_ ;

	return os;
}
/*
void DLFermion::create_Lz_basis(size_t k)
{
	if(k1==k)
		return;
	basis_u.create(N_phi, N_u, num_thread);
	//cout << basis_u << endl;
	basis_d.create(N_phi, N_d, num_thread);
	H_space.create_k1(basis_u, basis_d, k);
	H_space.hash();
	k1 = k;
	k2 = 0xFFFF;
	return;
}*/
/*

DLFermion A = DLFermion(8, 2, 2, 1);

extern "C" void init(size_t N_phi, size_t N_up, size_t N_down, int num_thread)
{A.reset(N_phi,N_up,N_down,num_thread);}

extern "C" void create_k1_space(size_t k1){A.create_k1_basis(k1);}
extern "C" void create_k2_space(size_t k2){A.create_k2_basis(k2);}

extern "C" size_t get_n_phi(){return A.get_N_phi();}
extern "C" size_t get_n_up(){return A.get_N_up();}
extern "C" size_t get_n_down(){return A.get_N_down();}
extern "C" size_t get_k1_dim(){return A.get_k1_dim();}
extern "C" size_t get_k2_dim(){return A.get_k2_dim();}

extern "C" size_t get_k2_n_data(unsigned int* indptr){return A.count_2symm(indptr);}

extern "C" size_t hamiltonian_generate(complex<double>* A1_list, complex<double>* A2_list, \
							complex<double>* A12_list, unsigned int* indptr, \
							unsigned int* indices, complex<double>* data)
{return A.generate_2symm(A1_list, A2_list, A12_list, indptr, indices, data);}

extern "C" void k2_to_k1_transformation(complex<double>* ket_k2, complex<double>* ket_k1)
{A.k2_to_k1_transformation(ket_k2, ket_k1);}

extern "C" void copy_basis(unsigned char *out, size_t i){A.copy_basis(out, i);}

extern "C" void occupation_number(complex<double>* ket, double *out)
{A.occupation_number(ket, out);}

extern "C" size_t to_3layers_trans(complex<double>* ket_k1, complex<double> *out)
{return A.to_3layers_trans(ket_k1, out);}*/