// Bipartition

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include "Basis.h"
#include "Bipartition.h"
#include "../time/time.h"
/*
long svd(Data_dtype* matrix, double* s, int m, int n)
{
	char jobu = 'N';
	char jobv = 'N';
	Data_dtype *u = NULL;
	Data_dtype *vt= NULL;
	int layout = LAPACK_ROW_MAJOR;
	Data_dtype *superb = (Data_dtype *)mkl_calloc(m*n, sizeof(Data_dtype), 64);
	double *s_t = s;

	auto info = LAPACKE_dgesvd(layout, jobu, jobv, m, n, matrix, n, s, u, 1, vt, 1, superb);
	if(info!=0)
		printf("Unconverged! Info=%lld\n", info);
	
	mkl_free(superb);
	return info;
}*/


//------------------------------------------------------------------------------
//----------------------------class  lzRealBasis--------------------------------
//------------------------------------------------------------------------------

lzRealBasis::lzRealBasis(int N_o, int N_e, int lz, const std::vector<Data_dtype>& alphalist)
{
	if(N_o<=1 || N_e<0 || N_e>2*N_o || lz<0)
	{
		auto res = boost::str( boost::format("(N_o = %d, N_e = %d, lz = %d)\n") % N_o % N_e % lz);
		throw std::runtime_error(res);
	}
	N_o_ = N_o;
	N_e_ = N_e;
	s_ = 0.5*N_o-0.5;
	lz_ = lz;
	N_tot_o_ = 2*N_o;

	int st = 0;
	size_t o = 0;
	st = fermionLzRealBasisCount(N_tot_o_, N_e_, lz_, N_e_, N_o_, o, 1.0, alphalist);
	if(st!=0)
	{
		auto res = boost::str( boost::format("fermionLzRealBasisCount\t%d\n") % st);
		throw std::runtime_error(res);
	}
	dim_ = o;
	auto temp = std::vector<Data_basis>(N_e_, 0);
	basis_ = std::vector<Data_basis>((size_t)dim_*N_e_, 0);
	coef_ = std::vector<Data_dtype>((size_t)dim_, 0.0);
	o = 0;
	st = fermionLzRealBasis(N_tot_o_, N_e_, lz_, temp.data(), basis_.data(), coef_.data(), N_e_, N_o_, o, 1.0, alphalist);
	if(st!=0||o!=dim_)
	{
		auto res = boost::str( boost::format("fermionLzRealBasis: (st=%d, o=%lu, dim=%lu)") % st % o % dim_);
		throw std::runtime_error(res);
	}
}
/*
std::ostream & operator<<(std::ostream & os, const lzRealBasis& t)
{
	int N_e = t.N_e();
	int N_o = t.N_o();
	auto it = t.basis().begin();
	bool printHashQ = (t.hash_.size()!=0);

	auto res = boost::format("N_o = %d, N_e = %d, lz = %d, dim = %lu ") % N_o % N_e % t.lz_ % t.dim_;
	os << "---------------------------------------------------------\n";
	os << res << (printHashQ?"HashTable:true":"HashTable:false") << std::endl;
	os << "---------------------------------------------------------\n";

	const lzBasis& fatherBasis = *static_cast<lzBasis const*>(&t);
	for(size_t i=0;i<t.dim_;++i)
	{
		os << fatherBasis.at(i);
		//printBasisLine(os, it, N_e, N_o, false);
		auto basenum = dict_num(it, N_e);
		if(printHashQ)
		{
			auto it = t.hash_.find(basenum);
			os << "\t" << basenum << "\t" << it->second;
		}
		else
		{
			os << "\t" << basenum << "\t" << i;
		}
		os << "\t" << boost::format("%.5e") % t.coef_.at(i) << std::endl;
		it += N_e;
	}
    //for(auto& it : t)
	//    std::cout << it;
	os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	os << res << (printHashQ?"HashTable:true":"HashTable:false") << std::endl;
	os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	return os;
}*/

// iterator

realBasisLine
lzRealBasis::operator[](size_t n)
{
	return realBasisLine(*this, n);
}

realBasisLine
lzRealBasis::at(size_t n) // check bound
{
	if( (n>=dim_)||(dim_==ULMAX_) )
	{
		auto res = boost::str( boost::format("dim=%lu, n=%lu\n") % dim_ % n);
		throw std::runtime_error(res);
	}
	return realBasisLine(*this,n);
}

realBasisLine lzRealBasis::cbegin() { return realBasisLine(*this,0); 	}
realBasisLine lzRealBasis::cend() 	{ return realBasisLine(*this,dim_); }
realBasisLine lzRealBasis::begin() 	{ return realBasisLine(*this,0); 	}
realBasisLine lzRealBasis::end() 	{ return realBasisLine(*this,dim_); }

std::ostream &  
lzRealBasis::printLine(size_t n) const // check bound
{
	//std::cout << "call realBasisLine\n";
	if( (n>=dim_)||(dim_==ULMAX_) )
	{
		auto res = boost::str( boost::format("dim=%lu, n=%lu\n") % dim_ % n);
		throw std::runtime_error(res);
	}
	return std::cout << realBasisLine(const_cast<lzRealBasis *>(this),n);
}

std::ostream & operator<<(std::ostream & os, const realBasisLine& t)
{
	//os << "operator<<(std::ostream & os, const realBasisLine& t)\n"; 
	os << static_cast<const basisLine&>(t);
	auto basis = *static_cast<lzRealBasis const*>(t.basis_ptr_);
	os << "\t" << basis.coef().at(t.it_);
	return os;
}
//---------------------------------------------------------------
//         basisOne
//---------------------------------------------------------------

basisOne::basisOne(const basisLine& obj) // copy obj
{
	auto lzbasis = obj.hilbert();
	int N_e = lzbasis.N_e();
	basis_ = std::vector<Data_basis>(N_e);
	for(auto i=0;i<N_e;++i)
		basis_[i] = obj[i];
	coef_ = obj.coef();
	base_num_ = dict_num( basis_.data(), N_e );
}

basisOne::basisOne(const basisLine& obj1, const basisLine& obj2)// obj1 + obj2;
{
	int N_e = obj1.N_e() + obj2.N_e();
	int N_temp = std::max(obj1[obj1.N_e()-1], obj2[obj2.N_e()-1])+1;
	auto temp = std::vector<Data_basis>(N_temp, 0);
	// insert: basis1->basis2
	for(auto i : obj2)
	{
		//std::cout << (signed)i << std::endl;
		temp[i] = 1;
	}
	//coef_ = static_cast<const realBasisLine*>(&obj1)->coef();
	int count = 0;
	int be = 0;
	int sign = 0;
	for(auto i : obj1)
	{
		if(temp[i])
			return;
		for(;be<=i;++be)
			count += temp[be];
		sign += count;
		temp[i] = 1;
		be = i+1;;
	}
	coef_ = obj1.coef() * obj2.coef();
	coef_ *= (sign&1)?-1.0:1.0;
	basis_ = std::vector<Data_basis>(N_e);
	//printf("N_o=%d\n", N_o);
	//for(auto i:temp)
	//	std::cout << (unsigned)i << " ";
	//std::cout << std::endl;
	binary_to_basis(temp.data(), basis_.data(), N_temp);
	base_num_ = dict_num( basis_.data(), N_e );
}

std::ostream & operator<<(std::ostream & os, const basisOne& obj)
{
	if(obj.basis().size()==0)
	{
		os << "NULL vector ";
	}
	else
	{
		for(auto i : obj.basis())
			os << (unsigned)i << " ";
		os << "\t" << obj.basenum() ;
	}
	os << boost::format("sign/coef: %+.5e\t")%obj.coef();// << std::endl;
	return os;
}

//---------------------------------------------------------------
//         basisOneMod
//---------------------------------------------------------------

std::ostream & operator<<(std::ostream & os, const basisOneMod& obj)
{
	if(obj.basis_.size()==0)
		os << "NULL vector ";
	else
	{
		printBasisLine(os, obj.basis().cbegin(), obj.basis().size(), obj.N_o(), false);
		os << "\t" << obj.basenum() << "\t";
		os << boost::format("sign/coef: %+.5e\t")%obj.coef();
	}
	//os << std::endl;
	return os;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int fermionLzRealBasis(Data_basis N_tot_o, Data_basis N_e, int Lz, Data_basis *temp, Data_basis *basis, Data_dtype *coef, const Data_basis N_p, const Data_basis N_o, size_t& o, Data_dtype initValue, const std::vector<Data_dtype>& alphalist)
//N_tot_o = 2*N_o
{
	if( std::abs(initValue)<1e-15 )
		return 0;
	if(N_e==0)
	{
		if(Lz!=0)
			return 0;
		Data_basis *basis_t = basis + (size_t)o * N_p;
		for(size_t j=0;j<N_p;++j)
		{
			basis_t[j] = temp[j];
			//printf("%d ",temp[j]);
		}
		coef[o] = initValue;
		//printf("\t%lu\n",dict_num(temp,N_p));
		++o;
		//getchar();
		return 0;
	}

	if(Lz<0)
		return 0;
	for(Data_basis i=N_e-1;i<N_tot_o;++i)
	{
		temp[(size_t)N_e - 1] = i;
		fermionLzRealBasis(i, N_e-1, Lz-(i%N_o), temp, basis, coef, N_p, N_o, o, initValue*alphalist[(i%N_o)], alphalist);
	}
	return 0;
}

int fermionLzRealBasisCount(int N_tot_o, int N_e, int lz, int N_p, int N_o, size_t& o, Data_dtype initValue, const std::vector<Data_dtype>& alphalist)
{
	if( std::abs(initValue)<1e-15 )
		return 0;
	if(N_e==0)
	{
		if(lz!=0)
			return 0;
		++o;
		return 0;
	}

	if(lz<0)
		return 0;
	for(Data_basis i=N_e-1;i<N_tot_o;++i)
		fermionLzRealBasisCount(i, N_e-1, lz-(i%N_o), N_p, N_o, o, initValue*alphalist[(i%N_o)], alphalist);
	return 0;
}

std::tuple<std::vector<Data_dtype>, std::vector<Data_dtype>> 
realCutACoef(int N_o, double theta)
{
	int N_phi = N_o - 1;
	//double theta = theta_p*M_PI*0.5;
	double cosp2 = std::pow( std::cos(0.5*theta) , 2);
	std::vector<Data_dtype> alphalist(N_o, 0.0);
	std::vector<Data_dtype> betalist(N_o, 0.0);
	for(int i=0;i<N_o;++i)
	{
		double temp = boost::math::ibeta(i+1, N_phi+1-i, cosp2);
		alphalist.at(i) = std::sqrt(temp);
		betalist.at(i) = std::sqrt(1-temp);
	}
	return std::make_tuple(alphalist, betalist);
}



//----------------------------------------------------------
//------------------------bipart----------------------------
//----------------------------------------------------------

std::tuple<int, int> bipart::lzRange(int Nl) const
{
	if(!checkNl(Nl))
	{
		auto res = boost::str( boost::format("(N_o = %d, N_e = %d, Nl = %d)\n") % No_ % Ne_ % Nl);
		throw std::runtime_error(res);
	}
	int lz_min = 0;
	int lz_max = 0;
	for(int i=Nl;i>1;i-=2)
	{
		lz_min += (Nl-i);
		lz_max += 2*No_ - 2 - (Nl-i);
	}
	if(Nl%2==1)
	{
		lz_min += int(Nl/2);
		lz_max += No_ - 1 - int(Nl/2);
	}
	if(Nl==0 || Nl==Ne_)
	{
		lz_min = 0;
		lz_max = 0;
	}
	return std::make_tuple(lz_min, lz_max);
}

bool bipart::checkNl(int Nl) const
{
	if(Nl<0 || Nl>Ne_)
		return false;
	return true;
}

reduceDensityMatrix bipart::reducedm(int Nl, int lzl, const std::vector<Data_dtype>& vec, const lzBasis& basis) const
{
	hu::timer hutm;
	if(Nl==0||Nl==basis.N_e())
		return this->reducedmBoundary(Nl, vec, basis);
	hutm.start("reducedm");
	auto [lz_min, lz_max] = this->lzRange(Nl);
	if(lzl<lz_min || lzl>lz_max)
	{
		auto res = boost::str( boost::format("(lz_min = %d, lz_max = %d, lzl = %d)\n") % lz_min % lz_max % lzl);
		throw std::runtime_error(res);
	}
	int L_z = basis.lz();
	hutm.start("lzRealBasis");
	auto basisl = lzRealBasis(No_, Nl, lzl, alphalist_);
	auto basisr = lzRealBasis(No_, Ne_-Nl, L_z-lzl, betalist_);
	hutm.stop("lzRealBasis");
	
	std::vector<Data_dtype> rdmData(basisl.dim()*basisr.dim(), 0.0);

	size_t dim_l = basisl.dim();
	size_t dim_r = basisr.dim();
	size_t UMAX = basis.hash_max();
	Data_dtype sum = 0.0;
	int numThreads = omp_get_max_threads();

	hutm.start("Combination");
	#pragma omp parallel for num_threads(numThreads) 
	for(int i=0;i<dim_l;++i)
	{
		for(int j=0;j<dim_r;++j)
		{
			auto refl = basisl[i];
			auto refr = basisr[j];
			auto bot = basisOne( refl, refr );
			//std::cout << bot%N_o;
			if(bot.nullVectorQ())
				continue;
			auto ind = basis.find(bot.basenum());
			//std::cout << refl << "+ ";
			//std::cout << refr << "= ";
			//std::cout << bot%No_ << "\t" << vec[ind] << "\t" << (bot.coef()*vec[ind])*(bot.coef()*vec[ind]) << std::endl;
			//sum += (bot.coef()*vec[ind])*(bot.coef()*vec[ind]);
			if(ind!=UMAX)
				rdmData[i*dim_r+j] = bot.coef()*vec[ind];
		}
	}
	hutm.stop("Combination");
	hutm.stop("reducedm");
	//printf("sum = %e\n", sum);
	return reduceDensityMatrix(dim_l, dim_r, rdmData);
}

reduceDensityMatrix bipart::reducedmBoundary(int Nl, const std::vector<Data_dtype>& vec, const lzBasis& basis) const
{
	hu::timer hutm;
	if(Nl!=0&&Nl!=basis.N_e())
	{
		auto res = boost::str( boost::format("In function reducedmBoundary: N_e=%d, Nl=%d\n")%basis.N_e() % Nl );
		throw std::runtime_error(res);
	}
	auto [lz_min, lz_max] = this->lzRange(Nl);
	int L_z = basis.lz();
	const std::vector<Data_dtype>& valuelist = (Nl==0)?betalist_:alphalist_;
	
	hutm.start("reducedmBoundary");
	auto basislr = lzRealBasis(No_, Ne_, L_z, valuelist);

	std::vector<Data_dtype> rdmData(basislr.dim(), 0.0);

	size_t dim = basislr.dim();
	for(int i=0;i<dim;++i)
	{
		auto ind = basis.find(basislr[i].basenum());
		rdmData[i] = basislr[i].coef()*vec[ind];
	}
	//printf("sum = %e\n", sum);
	hutm.stop("reducedmBoundary");
	if(Nl==0)
		return reduceDensityMatrix(1, dim, rdmData);
	return reduceDensityMatrix(dim, 1, rdmData);
}

//----------------------------------------------------------
//-----------------reduceDensityMatrix----------------------
//----------------------------------------------------------

std::vector<double> reduceDensityMatrix::singularValues(std::string method) const
{
	if(!initQ_)
		return std::vector<double>(0);
	bool wrongQ = false;

	hu::timer hutm;
	auto s = std::vector<double>(std::min(m_, n_), 0.0);
	auto matrix = rdm_;
	hutm.start(method);
	if( method=="gesvd" || method=="_gesvd" )
		mkl::dgesvd(matrix.data(), s.data(), m_, n_);
	else if( method=="gesdd" || method=="_gesdd" )
		mkl::dgesdd(matrix.data(), s.data(), m_, n_);
	else
		wrongQ = true;
	hutm.stop(method);

	if(wrongQ)
		throw std::runtime_error("Unrecognized method: " + method);
	return s;
}


std::vector<double> reduceDensityMatrix::eigenValues(std::string method) const
{
	if(!initQ_)
		return std::vector<double>(0);
	bool wrongQ = false;

	hu::timer hutm;
	long n = std::min(m_, n_);
	std::vector<Data_dtype> rdm(n*n, 0.0);
	const auto matrixA = matrixWrapper(const_cast<Data_dtype*>(rdm_.data()), m_, n_);
	auto matrixC = matrixWrapper(rdm.data(), n, n);
	bool transA = m_>n_;
	bool transB = !transA;
	
	hutm.start(method);
	hutm.start("dgemm");
	mkl::matrixMult(matrixA, matrixA, matrixC, transA, transB);
	hutm.stop("dgemm");

	auto w = std::vector<double>(std::min(m_, n_), 0.0);
	if( method=="dsyev" )
		mkl::dsyev (matrixC, w.data(), false, true);
	else if( method=="dsyevd" )
		mkl::dsyevd(matrixC, w.data(), false, true);
	else
		wrongQ = true;
	hutm.stop(method);

	if(wrongQ)
		throw std::runtime_error("Unrecognized method: " + method);
	return w;
}


std::vector<double> reduceDensityMatrix::diagonalLine(direction dir) const
{
	if(!initQ_)
		return std::vector<double>(0);
	size_t length = (dir==direction::ROW)?m_:n_;
	std::vector<double> resu( length, 0.0 );

	if(dir==direction::ROW)
	{
		for(size_t i=0;i<m_;++i)
			for(size_t j=0;j<n_;++j)
				resu[i] += rdm_[i*n_+j]*rdm_[i*n_+j];
	}
	else
	{
		for(size_t i=0;i<m_;++i)
			for(size_t j=0;j<n_;++j)
				resu[j] += rdm_[i*n_+j]*rdm_[i*n_+j];
	}
	return resu;
}


std::ostream& operator<<(std::ostream& os, const reduceDensityMatrix& obj)
{	
	size_t m = obj.dimRow();
	size_t n = obj.dimCol();
	const auto& rdm = obj.to_vec();
	auto res = boost::format("reduceDensityMatrix: dimRow = %lu, dimCol = %lu") % m % n;
	os << "---------------------------------------------------------\n";
	os << res << std::endl;
	os << "---------------------------------------------------------\n";

	for(size_t i=0;i<m;++i)
	{
		for(size_t j=0;j<n;++j)
		{
			std::cout << rdm[i*n+j] << "\t";
		}
		std::cout << std::endl;
	}
	return os;
}

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
res_ED(std::vector<Data_dtype> vec, std::vector<Data_dtype> & alphalist, std::vector<Data_dtype> & betalist, const lzBasis& basis, std::string method)
{
	hu::timer hutm;
	int N_e = basis.N_e();
	int N_o = basis.N_o();

	auto bp = bipart(basis, alphalist, betalist);

	double SA = 0;
	double norm_SA = 0;
	for(int Nl=0;Nl<=N_o;++Nl)
	{
		auto [lz_min, lz_max] = bp.lzRange(Nl);
		lz_max = std::min(lz_max, basis.lz());
		//std::cout << boost::format("lz_min = %d, lz_max = %d\n")%lz_min%lz_max;
		for(int lzl=lz_min;lzl<=lz_max;++lzl)
		{
			hutm.start("reducedm");
			auto rdm = bp.reducedm(Nl, lzl, vec, basis);
			hutm.stop("reducedm");
			hutm.start("entropy");
			auto [num_2, entropy] = rdm.entropy(method);
			hutm.stop("entropy");
			norm_SA += num_2;
			SA += entropy;
		}
		printf("Nl=%d, norm_SA = %lf, entropy = %lf\n", Nl, norm_SA, SA);
	}
	SA = SA/norm_SA+std::log(norm_SA);
	printf("Final EE: S_E/norm+log(norm) = %lf\n", SA);
	return std::make_tuple(SA, 1-norm_SA);
}

std::tuple<double, double, double, double> 
mutual_ED(double theta_p, std::vector<Data_dtype> vec, const lzBasis& basis, std::string method)
{
	hu::timer hutm;
	int N_e = basis.N_e();
	int N_o = basis.N_o();
	int N_phi = N_o - 1;
	std::vector<double> alphalist(N_o);
	std::vector<double> betalist(N_o);
	double theta = theta_p*M_PI*0.5;
	double cosp2 = std::pow( std::cos(theta) , 2);
	for(int i=0;i<N_o;++i)
	{
		double temp = boost::math::ibeta(i+1, N_phi+1-i, cosp2);
		betalist.at(i) = std::sqrt(temp);
		alphalist.at(i) = std::sqrt(1-temp);
	}

	auto [SA, errorA] = res_ED(vec, alphalist, betalist, basis, method);

	printf("S(A+B)\n");
	theta = (1-theta_p)*M_PI*0.5;
	double cosp2c = std::pow( std::cos(theta) , 2);
	for(int i=0;i<N_o;++i)
	{
		double temp = 1-boost::math::ibeta(i+1, N_phi+1-i, cosp2);
		temp += boost::math::ibeta(i+1, N_phi+1-i, cosp2c);
		//printfln("i=%d+1:\t%lf", i, temp);
		betalist.at(i) = std::sqrt(temp);
		alphalist.at(i) = std::sqrt(1-temp);
	}
	
	auto [SAB, errorAB] = res_ED(vec, alphalist, betalist, basis, method);

	return std::make_tuple(SA, errorA, SAB, errorAB);
}
