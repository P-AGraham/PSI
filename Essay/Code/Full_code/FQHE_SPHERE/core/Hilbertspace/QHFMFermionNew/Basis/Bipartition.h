

#pragma once

#include "Basis.h"
#include <boost/math/special_functions/beta.hpp>
#include <boost/format.hpp>
#include <tuple>
#include "../Algorithrm/gemm.h"
#include "../Algorithrm/matrixDecomp.h"
#include "../Algorithrm/printVector.h"
#include "../time/time.h"


class realBasisLine;
class lzRealBasis: public lzBasis
{
private:
    double threshHold_;
	std::vector<Data_dtype> coef_;
public:
	friend class realBasisLine;
	lzRealBasis() = default;
	lzRealBasis(int N_o, int N_e, int lz, const std::vector<Data_dtype>& alphaList);
	~lzRealBasis() = default;

	// 父类继承
	// value
	/*int N_o() const {return N_o_;};
	int N_e() const {return N_e_;};
	int N_phi() const {return N_phi_;};
	int N_tot_o() const {return N_tot_o_;};
	double s() const {return s_;};
	int lz() const {return lz_;};
	size_t dim() const {return dim_;};
	const std::vector<Data_basis>& basis() const {return basis_;};*/
	const std::vector<Data_dtype>& coef() const {return coef_;};

	// 父类继承
	// hash method
    /*const std::unordered_map<size_t, size_t>& generateIndex();
	const std::unordered_map<size_t, size_t>& hashTable() const {return hash_;};
	size_t find(size_t basenum) const;
	//size_t find(const_iterater it) const;
	size_t hash_max() const {return ULMAX_; };*/

	// iterator
	realBasisLine operator[](size_t n);
	realBasisLine at(size_t n); // check bound
	realBasisLine cbegin();
	realBasisLine cend();
	realBasisLine begin();
	realBasisLine end();
	virtual std::ostream & printLine(size_t n) const; // check bound
};

class realBasisLine: public basisLine
{
public:
	realBasisLine() = delete;
	realBasisLine(const lzBasis& basis, size_t it): basisLine(basis, it){};
	realBasisLine(const lzBasis* basis, size_t it): basisLine(basis, it){};

	virtual Data_dtype coef() const
	{	
		auto basis_ptr = static_cast<const lzRealBasis*>(basis_ptr_);
		return basis_ptr->coef()[it_];
	};

	friend std::ostream & operator<<(std::ostream & os, const realBasisLine& t);
};


class basisOne
{
protected:
	std::vector<Data_basis> basis_;
	size_t base_num_ = std::numeric_limits<size_t>::max();
	Data_dtype coef_ = 0;
public:
	basisOne() = default;
	basisOne(const basisLine& obj); // copy obj
	basisOne(const basisLine& obj1, const basisLine& obj2);
	// obj1 + obj2;
	~basisOne() = default;
	
	const std::vector<Data_basis>& basis() const {return basis_;};
	size_t basenum() const {return base_num_;};
	bool nullVectorQ() const {return basis_.size()==0;};
	Data_dtype coef() const {return coef_;};
};
std::ostream & operator<<(std::ostream & os, const basisOne& t);

// for print magic usage: basisOne%N_o;
class basisOneMod: public basisOne
{
private:
	int N_o_;
public:
	basisOneMod() = default;
	basisOneMod(const basisLine& obj, int N_o): basisOne(obj)
	{N_o_=N_o;}; // copy obj
	basisOneMod(const basisLine& obj1, const basisLine& obj2, int N_o): basisOne(obj1, obj2)
	{N_o_=N_o;};
	// obj1 + obj2;
	basisOneMod(const basisOne& obj, int N_o): basisOne(obj)
	{N_o_=N_o;}; // copy obj
	~basisOneMod() = default;
	
	int N_o() const {return N_o_;};
	friend std::ostream & operator<<(std::ostream & os, const basisOneMod& obj);
};

basisOneMod operator%(const basisOne& obj, int N_o)
{
	return basisOneMod(obj, N_o);
}
// modified magic usage: basisOne%N_o;
/*std::string operator%(const basisOne& obj, int N_o)
{
	return basisOneMod(obj, N_o);
}*/

int fermionLzRealBasisCount(int N_tot_o, int N_e, int lz, int N_p, int N_o, size_t& o, Data_dtype initValue, const std::vector<Data_dtype>& alphalist);

int fermionLzRealBasis(Data_basis N_tot_o, Data_basis N_e, int Lz, Data_basis *temp, Data_basis *basis, Data_dtype *coef, const Data_basis N_p, const Data_basis N_o, size_t& o, Data_dtype initValue, const std::vector<Data_dtype>& alphalist);
//N_tot_o = 2*N_o

std::tuple<std::vector<Data_dtype>, std::vector<Data_dtype>> 
realCutACoef(int N_o, double theta);











//----------------------------------------------------------
//------------------------bipart----------------------------
//----------------------------------------------------------
class reduceDensityMatrix;
class bipart
{
private:
	int No_;
	int No_tot_; // 2*No_
	int Ne_; // total number of electrons
	std::vector<Data_dtype> alphalist_;
	std::vector<Data_dtype> betalist_;

	reduceDensityMatrix reducedmBoundary(int Nl, const std::vector<Data_dtype>& vec, const lzBasis& basis) const;
	bool checkNl(int Nl) const;
public:
	bipart() = delete;
	bipart(const lzBasis& basis, const std::vector<Data_dtype>& alphalist, const std::vector<Data_dtype>& betalist)
	{
		No_ = basis.N_o();
		Ne_ = basis.N_e();
		No_tot_ = 2*No_;
		alphalist_ = alphalist;
		betalist_ = betalist;
	}
	~bipart() = default;
	//void generate_index(const Data_basis* basis)
	//template<class RandomAccessIterator>
	//std::tuple<double, double> entropy(Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist, int lz, size_t dim) const;
	
	//template<class RandomAccessIterator>
	reduceDensityMatrix reducedm(int Nl, int lzl, const std::vector<Data_dtype>& vec, const lzBasis& basis) const;

	std::tuple<int, int> lzRange(int Nl) const;
};

class reduceDensityMatrix // reduced density matrix
{
protected:
	size_t m_ = 0;
	size_t n_ = 0;
	std::vector<Data_dtype> rdm_;
	bool initQ_ = false;
	std::string checkMothod(std::string method) const
	{
		if(method=="_gesvd" || method=="gesvd")
			return "svd";
		if(method=="_gesdd" || method=="gesdd")
			return "svd";
		if(method=="dsyevd")
			return "eigen";
		if(method=="dsyev")
			return "eigen";
		return "";
	}
public:
	reduceDensityMatrix() = default;
	reduceDensityMatrix(size_t m, size_t n, std::vector<Data_dtype>& rdm_g)
	{
		m_ = m;
		n_ = n;
		rdm_ = rdm_g;
		initQ_ = true;
	}
	~reduceDensityMatrix() = default;

	const std::vector<Data_dtype>& to_vec() const {return rdm_;};
	size_t dimRow() const {return m_;};
	size_t dimCol() const {return n_;};
	size_t dimColumn() const {return n_;};

	std::vector<double> singularValues(std::string method="_gesvd") const;
	std::vector<double> eigenValues(std::string method="dsyev") const;
	std::tuple<double, double> entropy(std::string method="_gesvd") const
	{
		if(!initQ_)
			return std::make_tuple(0, 0);
		std::vector<double> s;
		if(checkMothod(method)=="svd" )
			s = this->singularValues(method);
		else if(checkMothod(method)=="eigen" )
			s = this->eigenValues(method);
		else
			throw std::runtime_error("reduceDensityMatrix::entropy: Unsupportted method: " + method);
		double norm_2 = 0.0;
		double entropy = 0.0;
		if(checkMothod(method)=="svd" )
		{
			for(auto i : s)
			{
				//printf("%lf\n", i);
				if(i>1e-20)
					entropy -= 2.0*i*i*std::log(i);
				norm_2 += i*i;
			}
		}
		else if(checkMothod(method)=="eigen" )
		{
			for(auto i : s)
			{
				if(i>1e-20)
					entropy -= i*std::log(i);
				norm_2 += i;
			}
		}
		return std::make_tuple(norm_2, entropy);
	}

	enum class direction{ROW, COLUMN};
	// the norm2 of every ROW of COLUMN
	std::vector<double> diagonalLine(direction dir) const;
};

std::ostream& operator<<(std::ostream& os, const reduceDensityMatrix& obj);


//----------------------------------------------------------
//-----------------Integrated RES method--------------------
//----------------------------------------------------------


/*
method
		if(method=="_gesvd" || method=="gesvd")
			return "svd";
		if(method=="_gesdd" || method=="gesdd")
			return "svd";
		if(method=="dsyevd")
			return "eigen";
		if(method=="dsyev")
			return "eigen";
*/
std::tuple<double, double> 
res_ED(std::vector<Data_dtype> vec, std::vector<Data_dtype> & alphalist, std::vector<Data_dtype> & betalist, const lzBasis& basis, std::string method);

std::tuple<double, double, double, double> 
mutual_ED(double theta_p, std::vector<Data_dtype> vec, const lzBasis& basis, std::string method);