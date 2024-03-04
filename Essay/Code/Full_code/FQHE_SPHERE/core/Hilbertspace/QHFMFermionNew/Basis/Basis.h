// Fermion Halperin State Basis.h

#ifndef __BASIS_H_
#define __BASIS_H_

#include "../QHFM.h"
#include <boost/format.hpp>

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<Lz_basis.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

size_t Binomial[50][50] = {0};

class basisLine;
class lzBasis
{
protected:
	int N_o_ = 1;
	int N_tot_o_ = 1;
	int N_e_ = 1;
	int N_phi_ = 1;
	double s_ = 0;
	int lz_;
	size_t dim_ = 0;
	std::vector<Data_basis> basis_;
	std::unordered_map<size_t, size_t> hash_;
	using const_iterater = std::vector<Data_basis>::const_iterator;
	const size_t ULMAX_ = std::numeric_limits<size_t>::max();
public:
	friend class basisLine;
	lzBasis() = default;
	lzBasis(int N_o, int N_e, int lz);
	~lzBasis() = default;

	// value
	int N_o() const {return N_o_;};
	int N_e() const {return N_e_;};
	int N_phi() const {return N_phi_;};
	int N_tot_o() const {return N_tot_o_;};
	double s() const {return s_;};
	int lz() const {return lz_;};
	size_t dim() const {return dim_;};
	const std::vector<Data_basis>& basis() const {return basis_;};

	// hash method
	const std::unordered_map<size_t, size_t>& generateIndex();
	const std::unordered_map<size_t, size_t>& hashTable() const {return hash_;};
	size_t find(size_t basenum) const;
	//size_t find(const_iterater it) const;
	size_t hash_max() const {return ULMAX_; };

	// iterator
	basisLine operator[](size_t n);
	basisLine at(size_t n); // check bound
	basisLine cbegin();
	basisLine cend();
	basisLine begin();
	basisLine end();
	virtual std::ostream & printLine(size_t n) const; // check bound
};
std::ostream & operator<<(std::ostream & os, const lzBasis& t);


// 不需要依赖友元类实现
// reference or const_iterator of single basis
class basisLine
{
protected:
	const lzBasis* basis_ptr_;
	size_t it_;
public:
	basisLine() = delete;
	basisLine(const lzBasis& basis, size_t it)
	{basis_ptr_=&basis;it_=it;};
	basisLine(const lzBasis* basis, size_t it)
	{basis_ptr_=basis;it_=it;};
	// assignment constructors are invalid
	basisLine(basisLine& obj) = default;
	basisLine(basisLine&& obj) = default;
	basisLine& operator=(const basisLine& obj) = default;
	basisLine& operator=(basisLine&& obj) = default;
	// copy method
	basisLine copy(){return basisLine(basis_ptr_, it_);}
	~basisLine() {it_=std::numeric_limits<size_t>::max();};
	
	// Hilbert space
	const lzBasis& hilbert() const {return *basis_ptr_;};
	int N_e() const {return basis_ptr_->N_e();};
	int N_o() const {return basis_ptr_->N_o();};
	int N_tot_o() const {return basis_ptr_->N_tot_o();};
	
	size_t basenum() const;
	virtual Data_dtype coef() const {return 1.0;};

	// == !=
	bool operator==(const basisLine& obj) const {
		return (basis_ptr_==obj.basis_ptr_)&&(it_==obj.it_);};
	bool operator!=(const basisLine& obj) const {
		return (basis_ptr_!=obj.basis_ptr_)||(it_!=obj.it_);};

	// ++ --
	basisLine& operator++() {++it_; return *this;};		// ++a
	basisLine operator++(int) {
		auto tmp = this->copy();
		++it_; return tmp;}; 	// a++
	basisLine& operator--() {--it_; return *this;};
	basisLine operator--(int) {
		auto tmp = this->copy();
		--it_; return tmp;};

	// + -
	basisLine& operator+(size_t n){it_ += n; return *this;};
	basisLine& operator-(size_t n){it_ -= n; return *this;};

	// [] at
	Data_basis operator[](size_t n) const 
	{ return basis_ptr_->basis_[it_*(basis_ptr_->N_e_)+n]; };
	Data_basis at(size_t n) const 
	{ return basis_ptr_->basis_.at(it_*(basis_ptr_->N_e_)+n); };
	const basisLine& operator*() const { return *this; };

	// iterator
	decltype(auto) cbegin() const {return basis_ptr_->basis_.cbegin()+it_*(basis_ptr_->N_e_);};
	decltype(auto) cend() const {return basis_ptr_->basis_.cbegin()+(it_+1)*(basis_ptr_->N_e_);};
	decltype(auto) begin() const {return basis_ptr_->basis_.cbegin()+it_*(basis_ptr_->N_e_);};
	decltype(auto) end() const {return basis_ptr_->basis_.cbegin()+(it_+1)*(basis_ptr_->N_e_);};
	
	friend std::ostream & operator<<(std::ostream & os, const basisLine& t);
};


template<class RandomIterator>
std::ostream & 
printBasisLine(std::ostream & os, const RandomIterator& it, int N_e, int N_o, bool ln_Q);
template<class RandomIterator>
std::ostream & 
printBasisLine(std::ostream & os, const RandomIterator& it, int N_e, int N_o);

int fermionLzBasis(Data_basis N_tot_o, Data_basis N_e, int Lz, Data_basis *temp,
					Data_basis *basis, const Data_basis N_p, const Data_basis N_o, size_t *o);//N_tot_o = 2*N_o;

int fermionLzBasisCount(Data_basis N_tot_o, Data_basis N_e, int Lz, const Data_basis N_p, const Data_basis N_o, size_t *o);

size_t binomial(size_t m, size_t n);//C_m^n;

template<class RandomIterator>
size_t dict_num(RandomIterator basis, size_t N_e);

size_t dict_num_binary(const Data_basis *basis, size_t N_e, size_t N_o_tot);

size_t gcd(size_t x, size_t y);


void basis_to_num(Data_basis *basis, size_t dim, size_t N_e, std::unordered_map<size_t, size_t> & hashT);

void binary_to_basis(Data_basis *ket_binary, Data_basis *ket_basis, size_t N_orbit);

void basis_to_binary(const Data_basis *ket_basis, Data_basis *ket_binary, size_t N_orbit, size_t N_e);

void print_basis(const Data_basis *basis, size_t N_e, size_t N_o);

void println_basis(const Data_basis *basis, size_t N_e, size_t N_o);

void print_basis_mod(const Data_basis *basis, size_t N_e, size_t N_o);

void println_basis_mod(const Data_basis *basis, size_t N_e, size_t N_o);

//-----------------------------------------------------------------------------------------------------------
//-------------------------<read.c>-------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

template<class T>
std::vector<T> read_vec(std::string &path);

#endif