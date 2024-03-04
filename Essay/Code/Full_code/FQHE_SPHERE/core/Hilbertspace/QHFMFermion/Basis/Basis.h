// Fermion Halperin State Basis.h

#ifndef __BASIS_H_
#define __BASIS_H_

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<Lz_basis.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

size_t Binomial[50][50] = {0};
//complex double exp_cc[50][50] = {0.0+0.0*I};


int fermionLzBasisCount(Data_basis N_tot_o, Data_basis N_e, int Lz, const Data_basis N_p, const Data_basis N_o, size_t *o);//N_tot_o = 2*N_o

int fermionLzBasis(Data_basis N_tot_o, Data_basis N_e, int Lz, Data_basis *temp,\
					Data_basis *basis, const Data_basis N_p, const Data_basis N_o, size_t *o);//N_tot_o = 2*N_o

size_t binomial(size_t m, size_t n);//C_m^n

size_t dict_num(const Data_basis *basis, size_t N_e);

size_t dict_num_binary(const Data_basis *basis, size_t N_e, size_t N_o_tot);

size_t gcd(size_t x, size_t y);

//int density(complex double *ket, unsigned char *k1_basis, double *out, \
			size_t N_e, size_t k1_dim);

//size_t *basis_to_num(Data_basis *basis, size_t dim, size_t N_e, MKL_INT *index);

void basis_to_num(Data_basis *basis, size_t dim, size_t N_e, std::unordered_map<size_t, size_t> & hashT);


void binary_to_basis(const Data_basis *ket_binary, Data_basis *ket_basis, size_t N_orbit); //
void basis_to_binary(const Data_basis *ket_basis, Data_basis *ket_binary, size_t N_orbit); //

void print_basis(const Data_basis *basis, size_t N_e, size_t N_o);

void println_basis(const Data_basis *basis, size_t N_e, size_t N_o);

bool projectionTest(const Data_basis *basis, size_t N_e, size_t N_o);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<Lz_Z2_basis.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

size_t base_translation(const Data_basis *basis, size_t offset, size_t N_orbit, size_t N_e, char *phase);

int fermionLzZ2BasisCount(Data_basis N_tot_o, Data_basis N_e, int Lz, const int Z2, Data_basis *temp,\
                      const Data_basis N_p, const Data_basis N_o, bool projectionQ, size_t *o);//N_tot_o = 2*N_o

int fermionLzZ2Basis(Data_basis N_tot_o, Data_basis N_e, int Lz, const int Z2, Data_basis *temp,\
                    Data_basis *basis, char* phase, const Data_basis N_p, const Data_basis N_o, bool projectionQ, size_t *o);

size_t find_Z2_index(const Data_basis *basis, const hash & hashT, size_t N_e, size_t N_o, int& n_find);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<Lz_Z2_PH_basis.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

void Z2_transformation(const Data_basis *basis, Data_basis *Z2basis, size_t N_o, size_t N_e, char *phase);

void PH_transformation(const Data_basis *basis, Data_basis *PHbasis, size_t N_o, size_t N_e, char *phase);

size_t find_Z2_index(const Data_basis *basis, const hash & hashT, const char* Z2_phase, size_t N_e, size_t N_o, int& length, char& phase);

size_t find_Z2PH_index(const Data_basis *basis, const hash & hashT, const hash & PHhashT, const char* Z2_phase,\
                    const char* PH_phase, const char* Z2PH_length, const size_t* PHindex, size_t N_e, size_t N_o, int& length, char& phase);

size_t find_Z2_index(const Data_basis *basis, const hash & hashT, const char* Z2_phase, size_t N_e, size_t N_o, int Z_2, int& length, char& phase);

size_t find_Z2PH_index(const Data_basis *basis, const hash & hashT, const hash & PHhashT, const char* Z2_phase,\
                    const char* PH_phase, const char* Z2PH_length, const size_t* PHindex, int Z_2, int PH, \
                        size_t N_e, size_t N_o, int& length, char& phase);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<density.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------


int density(Data_dtype *ket, Data_basis *basis, double *out, size_t N_e, size_t dim);

extern "C" int meausre_density(Data_dtype *ket, double *out);

int sigma_x(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, size_t N_e,  size_t N_o, size_t dim);

extern "C" int meausre_sigma_x(Data_dtype *ket, Data_dtype *out);


int Z2(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, size_t N_e, size_t dim);

extern "C" int meausre_Z2(Data_dtype *ket, Data_dtype *out);

extern "C" int meausre_PH(Data_dtype *ket, Data_dtype *out);

int PH(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, size_t N_e, size_t N_o, size_t dim);

extern "C" int meausre_PH_inZ2space(Data_dtype *ket, Data_dtype *out);

int PH(Data_dtype *ket, Data_basis *basis, Data_dtype *out, \
				const hash & hashT, char* Z2_phase, int Z_2, size_t N_e, size_t N_o, size_t dim);



//--------------------------------------------------------------------------------------------------------------------
//-------------------------<transformation.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" void LzZ2PH_to_Lz(const Data_dtype* LzZ2PH_ket, Data_dtype* Lz_ket, int num_thread);

void LzZ2PH_to_Lz_thread(const Data_dtype* LzZ2PH_ket, Data_dtype* Lz_ket, size_t row);


//-----------------------------------------------------------------------------------------------------------
//-------------------------<Bipartition.c>-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// Bipartition

size_t dict_num_offset(unsigned char *basis, size_t N_e, size_t offset);
/*
void orbitBipartition_thread(unsigned char *basis, struct bipartition *partition,\
						 	size_t N_e, size_t o_left);

void orbitBipartition(unsigned char *basis, struct bipartition *partition, size_t dim,\
					 	size_t N_e, size_t o_left, int num_thread);

void orbitBipartOrdering_Symmetry(unsigned char *basis, struct bipartition *partition, MKL_INT *rdmdim, \
									size_t occupy, size_t ky, size_t o_left);

void orbitBipartOrdering(unsigned char *basis, struct bipartition *partition, MKL_INT *rdmdim,\
						size_t o_left, int num_thread)//reduced density matrix;*/


int fermion1Basis(Data_basis N_tot_o, Data_basis N_e, std::vector<Data_basis>& temp,\
				std::vector<Data_basis>& basis, const Data_basis N_p, size_t *o);

struct bipart_inf
{
	size_t ind;
	//size_t ind;
	//Data_dtype coef = 0.0;
	//Data_dtype r_coef;
	//int lz;
};

class bipart
{
private:
	int No_;
	int No_tot_; // 2*No_
	int Ne_; // total number of electrons
	int Nl_; 
	int lz_min_ = 0;
	int lz_max_ = 0;
	int l_lz_;
	size_t mask_dim_;
	std::vector<Data_basis> l_mask_; // left occupation
	std::vector<char> sign_mask_;
	//std::vector<int> Lz_; 
	using lin_hash = std::unordered_map<size_t, bipart_inf>;
	std::vector<lin_hash> l_inf_; 
	std::vector<lin_hash> r_inf_; 
	std::unordered_map<size_t, std::tuple<size_t, size_t>> dim_;
	template<class RandomAccessIterator>
	void rdm_ij(size_t i, size_t j, int lz, size_t dimr, std::vector<Data_dtype>& rdm, Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist) const;

	template<class RandomAccessIterator>
	std::tuple<double, double> entropy(Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist, size_t dim) const;
public:
	bipart() = delete;
	template<class RandomAccessIterator>
	bipart(int N_o, int N_e, int N_l, int l_lz, RandomAccessIterator basis, size_t dim);
	~bipart() = default;
	//void generate_index(const Data_basis* basis)
	template<class RandomAccessIterator>
	std::tuple<double, double> entropy(Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist, int lz, size_t dim) const;
	int lz_min() const {return lz_min_;};
	int lz_max() const {return lz_max_;};
};

long svd(Data_dtype* matrix, double* s, int m, int n);

std::tuple<double, double> 
res_ED_ver1(int N_o, int N_e, std::vector<Data_dtype> vec, std::vector<Data_dtype> & alphalist, std::vector<Data_dtype> & betalist, const struct FQHE& A_);

std::tuple<double, double, double, double> 
mutual_ED_ver1(int N_o, int N_e, double theta_p, std::vector<Data_dtype> vec, const struct FQHE& A_);

#endif