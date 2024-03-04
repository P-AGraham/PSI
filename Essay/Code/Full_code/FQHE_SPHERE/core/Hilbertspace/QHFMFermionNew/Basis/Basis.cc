//  3Body_Fermion_k1_basis.h

#include "Basis.h"

//---------------------------------------------------------------
//         lzBasis
//---------------------------------------------------------------

lzBasis::
lzBasis(int N_o, int N_e, int lz)
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
	st = fermionLzBasisCount(N_tot_o_, N_e_, lz_, N_e_, N_o_, &o);
	if(st!=0)
	{
		auto res = boost::str( boost::format("fermionLzBasisCount\t%d\n") % st);
		throw std::runtime_error(res);
	}
	dim_ = o;
	auto temp = std::vector<Data_basis>(N_e_, 0);
	basis_ = std::vector<Data_basis>((size_t)o*N_e_, 0);
	o = 0;
	st = fermionLzBasis(N_tot_o_, N_e_, lz_, temp.data(), basis_.data(), N_e_, N_o_, &o);
	if(st!=0||o!=dim_)
	{
		auto res = boost::str( boost::format("fermionLzBasis: (st=%d, o=%lu, dim=%lu)") % st % o % dim_);
		throw std::runtime_error(res);
	}
}

const std::unordered_map<size_t, size_t>& lzBasis::
generateIndex()
{
	auto basis_t = basis_.data();
	for(size_t i=0;i<dim_;++i)
	{
		size_t basenum = dict_num(basis_t, N_e_);
		hash_.insert({basenum, i});
		basis_t += N_e_;
	}
	return hash_;
}

size_t lzBasis::find(size_t basenum) const
{
	auto it = hash_.find(basenum);
	if(it==hash_.end())
		return ULMAX_;
	return it->second;
}

std::ostream & operator<<(std::ostream & os, const lzBasis& t)
{
	int N_e = t.N_e();
	int N_o = t.N_o();
	auto it = t.basis().begin();
	bool printHashQ = (t.hashTable().size()!=0);

	auto res = boost::format("N_o = %d, N_e = %d, lz = %d, dim = %lu ") % N_o % N_e % t.lz() % t.dim();
	os << "---------------------------------------------------------\n";
	os << res << (printHashQ?"HashTable:true":"HashTable:false") << std::endl;
	os << "---------------------------------------------------------\n";

	for(size_t i=0;i<t.dim();++i)
	{
		t.printLine(i); 
		//printBasisLine(os, it, N_e, N_o, false);
		auto basenum = dict_num(it, N_e);
		if(printHashQ)
		{
			auto itf = t.hashTable().find(basenum);
			os << "\t" << basenum << "\t" << itf->second;
		}
		else
		{
			os << "\t" << basenum << "\t" << i;
		}
		os << std::endl;
		it += N_e;
	}
    //for(auto& it : t)
	//    std::cout << it;
	os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	os << res << (printHashQ?"HashTable:true":"HashTable:false") << std::endl;
	os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	return os;
}


// iterator

basisLine 
lzBasis::operator[](size_t n)
{
	return basisLine(this, n);
}

basisLine 
lzBasis::at(size_t n) // check bound
{
	if( (n>=dim_)||(dim_==ULMAX_) )
	{
		auto res = boost::str( boost::format("dim=%lu, n=%lu\n") % dim_ % n);
		throw std::runtime_error(res);
	}
	return basisLine(this,n);
}

basisLine lzBasis::cbegin() { return basisLine(this,0); 	}
basisLine lzBasis::cend() 	{ return basisLine(this,dim_); }
basisLine lzBasis::begin() 	{ return basisLine(this,0); 	}
basisLine lzBasis::end() 	{ return basisLine(this,dim_); }

std::ostream & 
lzBasis::printLine(size_t n) const// check bound
{
	//std::cout << "call basisLine\n";
	if( (n>=dim_)||(dim_==ULMAX_) )
	{
		auto res = boost::str( boost::format("dim=%lu, n=%lu\n") % dim_ % n);
		throw std::runtime_error(res);
	}
	return std::cout << basisLine( const_cast<lzBasis*>(this),n);
}


//---------------------------------------------------------------
//         basisLine
//---------------------------------------------------------------

size_t basisLine::basenum() const
{
	int N_e = basis_ptr_->N_e();
	return dict_num(basis_ptr_->basis().cbegin()+it_*N_e, N_e);
}

std::ostream & 
operator<<(std::ostream & os, const basisLine& t)
{
	//os << "operator<<(std::ostream & os, const basisLine& t)\n"; 
	int N_e = t.basis_ptr_->N_e();
	int N_o = t.basis_ptr_->N_o();
	auto it = t.basis_ptr_->basis().begin() + t.it_*N_e;
	//os << boost::format( "Line %lu: (Ne=%d, No=%d)" ) % t.it_ % N_e % N_o;
	return printBasisLine(os, it, N_e, N_o, false);
}

template<class RandomIterator>
std::ostream & 
printBasisLine(std::ostream & os, const RandomIterator& it, int N_e, int N_o, bool ln_Q)
{
	for(auto i=0;i<N_e;++i)
	{
		if( it[i]<N_o )
			os << (unsigned)it[i] << " ";
		else
			os << "\033[1;31m" << (unsigned)it[i]%N_o << "\033[0m" << " ";
	}
	if(ln_Q)
		os << std::endl;
	return os;
}
template<class RandomIterator>
std::ostream & 
printBasisLine(std::ostream & os, const RandomIterator& it, int N_e, int N_o)
{
	return printBasisLine(os, it, N_e, N_o, true);
}
















//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int fermionLzBasis(Data_basis N_tot_o, Data_basis N_e, int Lz, Data_basis *temp,\
					Data_basis *basis, const Data_basis N_p, const Data_basis N_o, size_t *o)//N_tot_o = 2*N_o
{
	if(N_e==0)
	{
		if(Lz!=0)
			return 0;
		Data_basis *basis_t = basis + (size_t)o[0] * N_p;
		for(size_t j=0;j<N_p;++j)
		{
			basis_t[j] = temp[j];
			//printf("%d ",temp[j]);
		}
		//printf("\t%lu\n",dict_num(temp,N_p));
		++o[0];
		//getchar();
		return 0;
	}

	if(Lz<0)
		return 0;
	for(Data_basis i=N_e-1;i<N_tot_o;++i)
	{
		temp[(size_t)N_e - 1] = i;
		fermionLzBasis(i, N_e-1, Lz-(i%N_o), temp, basis, N_p, N_o, o);
	}
	return 0;
}


int fermionLzBasisCount(Data_basis N_tot_o, Data_basis N_e, int Lz, const Data_basis N_p, const Data_basis N_o, size_t *o)
{
	if(N_e==0)
	{
		if(Lz!=0)
			return 0;
		++o[0];
		return 0;
	}

	if(Lz<0)
		return 0;
	for(Data_basis i=N_e-1;i<N_tot_o;++i)
		fermionLzBasisCount(i, N_e-1, Lz-(i%N_o), N_p, N_o, o);
	return 0;
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

template<class RandomIterator>
size_t dict_num(RandomIterator basis, size_t N_e)
{
	size_t num = 0;
	for(size_t i=0;i<N_e;++i) // particle number
		num += binomial(basis[i], i+1);
	return num;
}

size_t dict_num_binary(const Data_basis *basis, size_t N_e, size_t N_o_tot)
{
	size_t num = 0, i=0;
	for(size_t j=0;j<N_o_tot;++j)
	{
		if(basis[j])
		{
			++i;
			num += binomial(j, i);
			//printf("binomial(%lu, %lu) + ", j, i);
		}
	}
		//printf("\n");
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


void basis_to_num(Data_basis *basis, size_t dim, size_t N_e, std::unordered_map<size_t, size_t> & hashT)
{
	Data_basis *basis_t = basis;
	size_t base_num;

	for(size_t i=0;i<dim;++i)
	{
		base_num = dict_num(basis_t, N_e);
		hashT.insert({base_num, i});
		basis_t += (size_t)N_e;
	}
	return;
}

void binary_to_basis(Data_basis *ket_binary, Data_basis *ket_basis, size_t N_orbit)
{
    // translate binary -> occupy
    size_t k = 0;
    for(size_t i=0;i<N_orbit;++i)
    {
        if(ket_binary[i])
        {
            ket_basis[k] = i;
            ++k;
        }
    }
}

void basis_to_binary(const Data_basis *ket_basis, Data_basis *ket_binary, size_t N_orbit, size_t N_e)
{
    // translate occupy -> binary
    memset(ket_binary, 0, (size_t)N_orbit * sizeof(Data_basis));
    for(auto i=0;i<N_e;++i)
    	ket_binary[ket_basis[i]] = 1;
}

void print_basis(const Data_basis *basis, size_t N_e, size_t N_o)
{
	for(size_t j=0;j<N_e;++j)
	{
		if(basis[j]<N_o)
			printf("%d ",basis[j]);
		else
			printf("\033[1;31m%d \033[0m",basis[j]);
	}
}

void println_basis(const Data_basis *basis, size_t N_e, size_t N_o)
{
	print_basis(basis, N_e, N_o);
	printf("\n");
}

void print_basis_mod(const Data_basis *basis, size_t N_e, size_t N_o)
{
	for(size_t j=0;j<N_e;++j)
	{
		if(basis[j]<N_o)
			printf("%d ",basis[j]);
		else
			printf("\033[1;31m%lu \033[0m",basis[j]%N_o);
	}
}

void println_basis_mod(const Data_basis *basis, size_t N_e, size_t N_o)
{
	print_basis_mod(basis, N_e, N_o);
	printf("\n");
}