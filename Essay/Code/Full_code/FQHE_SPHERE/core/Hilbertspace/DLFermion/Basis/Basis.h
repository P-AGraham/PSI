
//--------------------------------------------------------------------------------------------------------------------
//-------------------------<tools.cpp>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

size_t Binomial[100][100] = {0};
extern "C" size_t binomial(size_t m, size_t n);

extern "C" int fermion_basis(size_t N_phi, size_t N_e, unsigned char *temp,\
					unsigned char *basis, const size_t N_p, size_t *o);

size_t gcd(size_t x, size_t y);

size_t dict_num(const unsigned char* basis, size_t N_e);

size_t dict_num_minus_phi(const unsigned char* basis, size_t N_e, size_t N_phi);

size_t base_translation(const unsigned char* basis, size_t offset,\
						size_t N_phi, size_t N_e, char *phase);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<k1_basis.cpp>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

//inline unsigned char cal_moment_thread(unsigned char *basis, unsigned char N_e, unsigned char N_phi);
//inline unsigned char cal_moment_thread(unsigned char *basis, unsigned char N_e, unsigned char N_phi, unsigned char N_m);
//--------------------------------------------------------------------------------------------------------------------
//-------------------------<basis_trans.cpp>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

//bool tlsyers_Q(const unsigned char *basis, size_t N_phi, size_t N_u, size_t N_d);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<class basis_array>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
// Lz_basis.cpp
class basis_array
{
private:
	unsigned char *basis;
	unsigned char *moment;
	size_t N_phi;
	size_t N_e;
	size_t dim;
public:
	basis_array(){basis=NULL;moment=NULL;N_e=0;N_phi=0;dim=0;}
	~basis_array();
	void reset();
	size_t create(size_t N_phi, size_t N_ele, int num_thread, size_t N_mod=0);
	const unsigned char *operator[](size_t i){return (i<dim)?(basis+i*N_e):NULL;}
	unsigned char get_k(size_t i){return (i<dim)?moment[i]:0xFF;}
	size_t get_N_e(){return N_e;}
	size_t get_N_phi(){return N_phi;}
	size_t get_dim(){return dim;}
	friend std::ostream & operator<<(std::ostream & os, const basis_array & b);
};

basis_array::~basis_array()
{
	if(basis!=NULL)
		delete [] basis;
	if(moment!=NULL)
		delete [] moment;
	N_e = 0;
	N_phi = 0;
	dim = 0;

	//std::cout << "destroy basis_array at " << this << std::endl;
}

void basis_array::reset()
{
	if(basis!=NULL)
		delete [] basis;
	basis = NULL;
	if(moment!=NULL)
		delete [] moment;
	moment = NULL;
	N_e = 0;
	N_phi = 0;
	dim = 0;
}
/*
//--------------------------------------------------------------------------------------------------------------------
//-------------------------<class basis>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

class basis
{
private:
	const unsigned char *u; //init
	const unsigned char *d; //init
public:
	size_t base_num;		//set_base_num(hash)
	size_t index;			//init
	size_t first = 0;
	char period =0 ;
	char sn = 0;
	char phase = 0;
	basis(){u=NULL;d=NULL;index=0;sn = 0;}
	void init(const unsigned char *up, const unsigned char *down, size_t ind){u=up; d=down;index=ind;}
	~basis(){};
	const unsigned char *up(){return u;}
	const unsigned char *down(){return d;}
	//void set_base_num(size_t t){base_num=t;}
	void print(size_t N_up, size_t N_down);
	//friend std::ostream & operator<<(std::ostream & os, const basis & b);
};

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<class basis_symm>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

class basis_symm
{
private:
	basis *k1_basis = NULL;
	basis *ket;
	size_t N_u;
	size_t N_d;
	size_t N_phi;
	size_t k1_dim;
	size_t k2_dim = 0;
	size_t k1 = 0xFFFF;
	size_t k2 = 0xFFFF;
	size_t stride;
	size_t p;
	size_t q;
	size_t pqgcd;
	size_t translational_symm_Q = 0;
	size_t *k2_basis = NULL;
	unordered_map<size_t, size_t> k1_hash;
	unordered_map<size_t, size_t> k2_hash;
public:
	basis_symm(){k1_basis=NULL;k2_basis = NULL;}
	~basis_symm();
	void reset();

	friend std::ostream & operator<<(std::ostream & os, basis_symm & b);
	void show(){cout << "k1 = " << k1 << "\tk2 = " << k2;\
				cout << "\tdim1 = " << k1_dim << "\tdim2 = " << k2_dim << endl;}

	size_t get_k1_dim(){return k1_dim;}
	size_t get_k2_dim(){return k2_dim;}

	void operator++(){++ket;}
	void operator--(){--ket;}
	void set(size_t i){ket=k1_basis+i;}
	size_t find(size_t key){return k1_hash[key];}
	//size_t find_k2(size_t key){return k2_hash[key];}
	size_t key(){ket->base_num=dict_num(ket->up(),N_u)*stride+dict_num(ket->down(),N_d);return ket->base_num;}
	size_t find(){return k1_hash[this->key()];} // cal the key of *ket and return the value k1_hash[key];
	//size_t find_k2(){return k2_hash[this->key()];}

	void create_k1(basis_array &u, basis_array &d, size_t k);
	void translational_symm();
	void create_k2(size_t k);
	void hash();

	void copy_basis(unsigned char *out, size_t i);
	void occupation_number(complex<double>* ket, double *out, int num_thread);

	size_t element_count(unsigned int* indptr, int num_thread);
	void element_count_thread(unsigned int row, unsigned int* indptr);
	size_t element_generate_thread(unsigned int row, \
			complex<double>* A1_list, complex<double>* A2_list, complex<double>* A12_list,\
			unsigned int* indptr, unsigned int* indices, complex<double>* data);
	size_t element_generate(complex<double>* A1_list, \
			complex<double>* A2_list, complex<double>* A12_list, unsigned int* indptr, \
			unsigned int* indices, complex<double>* data, int num_thread);

	void k2_to_k1_thread(size_t i, complex<double>* ket_k2, complex<double>* ket_k1);
	size_t to_3layers_trans(complex<double>* ket_k1, complex<double> *out, int num_thread);
};*/
