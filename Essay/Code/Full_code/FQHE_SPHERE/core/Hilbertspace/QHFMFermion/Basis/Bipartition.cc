// Bipartition

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include "Basis.h"
#include "../time/time.h"

size_t dict_num_offset(unsigned char *basis, size_t N_e, size_t offset)
{
	unsigned int num = 0;
	for(size_t i=0;i<N_e;++i) // particle number
		num += binomial(basis[i]-offset, i+1);
	return num;
}
/*
void orbitBipartition_thread(unsigned char *basis, struct bipartition *partition,\
						 	size_t N_e, size_t o_left)
{
	size_t p_left = 0;
	size_t moment = 0;
	size_t i = 0;
	while((basis[i]<o_left) && (i<N_e))
	{
		moment += basis[i]%A.N_o;
		++i;
	}
	partition->occupy = i;
	partition->Lz = moment;
	partition->left = dict_num(basis, i);
	partition->right = dict_num_offset(basis+i, N_e-i, o_left);
}

void orbitBipartition(unsigned char *basis, struct bipartition *partition, size_t dim,\
					 	size_t N_e, size_t o_left, int num_thread)
{
	# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<dim;++i)
		orbitBipartition_thread(basis+i*N_e, partition+i, N_e, o_left);
}

void orbitBipartOrdering_Symmetry(unsigned char *basis, struct bipartition *partition, MKL_INT *rdmdim, \
									size_t occupy, size_t ky, size_t o_left)
{
	MKL_INT *l_index = (MKL_INT *)malloc((size_t)binomial(o_left, occupy) * sizeof(MKL_INT));
	MKL_INT *r_index = (MKL_INT *)malloc((size_t)binomial(A.N_phi-o_left, A.N_e-occupy) * sizeof(MKL_INT));
	memset(l_index, -1, (size_t)binomial(o_left, occupy) * sizeof(MKL_INT));
	memset(r_index, -1, (size_t)binomial(A.N_phi-o_left, A.N_e-occupy) * sizeof(MKL_INT));
	MKL_INT left, right = 0;
	left = 0;
	for (size_t i=0;i<A.dim;++i)
	{
		if(partition[i].occupy!=occupy || partition[i].ky!=ky)
			continue;
		if(l_index[partition[i].left]==0xFFFFFFFF)
		{
			l_index[partition[i].left] = left;
			++left;
		}
		if(r_index[partition[i].right]==0xFFFFFFFF)
		{
			r_index[partition[i].right] = right;
			++right;
		}
		partition[i].left = l_index[partition[i].left];
		partition[i].right = r_index[partition[i].right];
	}
	free(l_index);
	free(r_index);
	rdmdim[(occupy*A.N_phi+ky)*2] = left;
	rdmdim[(occupy*A.N_phi+ky)*2+1] = right;
	//printf("%d\t%d:\t%u\t%u\n",occupy,ky,left,right);
}

void orbitBipartOrdering(unsigned char *basis, struct bipartition *partition, MKL_INT *rdmdim,\
						size_t o_left, int num_thread)//reduced density matrix
{
	# pragma omp parallel for num_threads(num_thread)
	for(size_t i=0;i<=o_left;++i)
		for (size_t ky=0;ky<A.N_phi;++ky)
			orbitBipartOrdering_Symmetry(basis, partition, rdmdim, i, ky, o_left);
}*/

int fermion1Basis(Data_basis N_tot_o, Data_basis N_e, std::vector<Data_basis>& temp,\
				std::vector<Data_basis>& basis, const Data_basis N_p, size_t *o)
{
	if(N_e==0)
	{
		auto basis_t = basis.begin() + (size_t)o[0] * N_p;
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

	for(Data_basis i=N_e-1;i<N_tot_o;++i)
	{
		temp[(size_t)N_e - 1] = i;
		fermion1Basis(i, N_e-1, temp, basis, N_p, o);
	}
	return 0;
}

#ifndef BIPART_DEBUG
//#define BIPART_DEBUG
#endif

template<class RandomAccessIterator>
bipart::
bipart(int N_o, int N_e, int N_l, int l_lz, RandomAccessIterator basis, size_t dim)
{
	hu::timer hutm;
	No_ = N_o;
	No_tot_ = 2*N_o;
	Ne_ = N_e;
	Nl_ = N_l;
	l_lz_ = l_lz;
	//printf("N_o=%d, N_e=%d, N_l=%d, l_lz=%d\n", N_o, N_e, N_l, l_lz);

	// lz range
	lz_min_ = 0;
	lz_max_ = 0;
	for(int i=Nl_;i>1;i-=2)
	{
		lz_min_ += (Nl_-i);
		lz_max_ += 2*No_ - 2 - (Nl_-i);
	}
	if(Nl_%2==1)
	{
		lz_min_ += int(Nl_/2);
		lz_max_ += No_ - 1 - int(Nl_/2);
	}
	if(Nl_==0 || Nl_==N_e)
	{
		lz_min_ = 0;
		lz_max_ = 0;
		return;
	}

	// generate mask
	mask_dim_ = binomial(Ne_, Nl_);
	l_mask_ = std::vector<Data_basis>(mask_dim_*(size_t)Nl_);
	std::vector<Data_basis> temp(Nl_);
	size_t o = 0;
	fermion1Basis(Ne_, Nl_, temp, l_mask_, Nl_, &o);
	#ifdef BIPART_DEBUG
	for(auto i=0;i<mask_dim_;++i)
	{
		for(auto j=0;j<Nl_;++j)
		{
			printf("%d ", l_mask_[i*Nl_+j]);
		}
		printf("\n");
	}
	#endif

	//printf("lz range :(%d, %d)\n", lz_min_, lz_max_);

	// bipartition
	std::vector<size_t> l_count((lz_max_-lz_min_+1), 0);
	std::vector<size_t> r_count((lz_max_-lz_min_+1), 0);
	l_inf_ = std::vector<lin_hash>(lz_max_-lz_min_+1);
	r_inf_ = std::vector<lin_hash>(lz_max_-lz_min_+1);

	hutm.start("sign");
	sign_mask_ = std::vector<char>(mask_dim_, 1.0);
	int numThreads = omp_get_max_threads();
	#pragma omp parallel for num_threads(numThreads)
	for(size_t i=0;i<mask_dim_;++i)
	{
		int sign_c = 0;
		for(size_t j=0;j<Nl_;++j)
		{
			sign_c += (l_mask_[i*Nl_+j]-j);
		}
		sign_mask_[i] = (sign_c&1==1)?(-1):1;
	}
	hutm.stop("sign");

	//printf("\n");
	for(size_t i=0;i<dim;++i)
	{
		#ifndef BIPART_DEBUG
		//printf("\r%lu/%lu", i, dim);
		#endif
		auto basis_b = basis + i*Ne_;
		auto basis_e = basis + (i+1)*Ne_;
		#ifdef BIPART_DEBUG
		println_basis_mod(basis_b, Ne_, No_);
		#endif
		for(size_t j=0;j<mask_dim_;++j)
		{
			int lz_l = 0;
			size_t l_key = 0;
			size_t r_key = 0;
			auto mask_b = l_mask_.cbegin() + j*Nl_;
			auto mask_e = l_mask_.cbegin() + (j+1)*Nl_;
			int l_c = 0;
			int r_c = 0;

			bipart_inf l_data;
			bipart_inf r_data;

			//hutm.start("count_lz");
			for(auto b_it=mask_b;b_it!=mask_e;b_it++)
			{
				lz_l += basis_b[*b_it]%No_;
			}
			//hutm.stop("count_lz");
			if(l_lz_>0 && lz_l!=l_lz_)
				continue;
			
			//hutm.start("make_key_value");
			for(auto b_it=basis_b;b_it!=basis_e;b_it++)
			{
				if(*b_it==basis_b[*mask_b] && mask_b!=mask_e)// left
				{
					//lz_l += (*b_it)%No_;
					l_key += binomial(*b_it, l_c+1);
					++l_c;
					++mask_b; 
				}
				else// right
				{
					r_key += binomial(*b_it, r_c+1);
					++r_c;
				}
			}
			//hutm.stop("make_key_value");
			// add to unordered_map

			//hutm.start("make_hash");
			#ifdef BIPART_DEBUG
			mask_b = l_mask_.cbegin() + j*Nl_;
			mask_e = l_mask_.cbegin() + (j+1)*Nl_;
			printf("\t left: ");
			for(auto b_it=mask_b;b_it!=mask_e;b_it++)
			{
				printf("%d ", basis_b[*b_it]);
			}
			#endif

			auto lit = l_inf_.at(lz_l-lz_min_).find(l_key);
			if(lit==l_inf_.at(lz_l-lz_min_).end())
			{
				l_data.ind = l_count.at(lz_l-lz_min_);
				l_count.at(lz_l-lz_min_) += 1;
				l_inf_.at(lz_l-lz_min_).insert({l_key, l_data});
				
				#ifdef BIPART_DEBUG
				printf("\t\033[4;31madded\033[0m, Lz=%d, dict_num=%lu, ind=%lu, coef=%e\n", lz_l, l_key, l_data.ind, l_data.coef);
				#endif
			}
			#ifdef BIPART_DEBUG
			else
				printf("\texist, Lz=%d, dict_num=%lu, ind=%lu, coef=%e\n", lz_l, lit->first, lit->second.ind, lit->second.coef);
			#endif

			#ifdef BIPART_DEBUG
			printf("\t right: ");
			for(auto b_it=basis_b;b_it!=basis_e;b_it++)
			{
				if(*b_it==basis_b[*mask_b] && mask_b!=mask_e)// left
					++mask_b;
				else// right
					printf("%d ", *b_it);
			}
			#endif

			auto rit = r_inf_.at(lz_l-lz_min_).find(r_key);
			if(rit==r_inf_.at(lz_l-lz_min_).end())
			{
				r_data.ind = r_count.at(lz_l-lz_min_);
				r_count.at(lz_l-lz_min_) += 1;
				r_inf_.at(lz_l-lz_min_).insert({r_key, r_data});
				
				#ifdef BIPART_DEBUG
				printf("\t\033[4;31madded\033[0m, Lz=%d, dict_num=%lu, ind=%lu, coef=%lf\n", lz_l, r_key, r_data.ind, r_data.coef);
				#endif
			}
			#ifdef BIPART_DEBUG
			else
				printf("\texist, Lz=%d, dict_num=%lu, ind=%lu, coef=%lf\n", lz_l, rit->first, rit->second.ind, rit->second.coef);
			#endif
			//hutm.stop("make_hash");
		}
	}
	#ifdef BIPART_DEBUG
	{
	int c = 0;
	for(const auto& i : l_inf_)
	{
		printf("left: L_z = %d\n", lz_min_+c);
		for(const auto j : i)
		{
			printf("\tdict_num=%lu, ind=%lu, coef=%lf\n", j.first, j.second.ind, j.second.coef);
		}
		++c;
	}c=0;
	for(const auto& i : r_inf_)
	{
		printf("right: L_z = %d\n", lz_min_+c);
		for(const auto& j : i)
		{
			printf("\tdict_num=%lu, ind=%lu, coef=%lf\n", j.first, j.second.ind, j.second.coef);
		}
		++c;
	}
	}
	#endif
	{
	size_t l_t = 0;
	size_t r_t = 0;
	for(auto i=0;i<l_count.size();++i )
	{
		l_t += l_count.at(i);
		r_t += r_count.at(i);
		//printf("%d:\t(%lu\t+\t%lu\t=%lu)\n",lz_min_+i ,l_count.at(i), r_count.at(i), l_count.at(i) + r_count.at(i));
		dim_.insert({i+lz_min_, std::make_tuple(l_count.at(i),r_count.at(i))});
	}//printf("l_dim = %lu\t r_dim = %lu\t total = %lu\n", l_t, r_t, dim);
	}
	//hutm.showAllSorted();
}

template<class RandomAccessIterator>
std::tuple<double, double> bipart::
entropy(Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist, int lz, size_t dim) const
{
	hu::timer hutm;
	if(Nl_==0 || Nl_==Ne_)
	{
		return entropy(vector, basis, alphalist, betalist, dim);
	}

	auto dit = dim_.find(lz);
	if(dit==dim_.end())
		return std::make_tuple(0.0, 0.0);
	auto diml = std::get<0>(dit->second);
	auto dimr = std::get<1>(dit->second);
	if(diml*dimr==0)
		return std::make_tuple(0.0, 0.0);
	//printf("lz = %d, diml = %lu, dimr = %lu\n", lz, diml, dimr);
	
	std::vector<Data_dtype> rdm(diml*dimr, 0.0);
	const lin_hash& lhash = l_inf_.at(lz-lz_min_);
	const lin_hash& rhash = r_inf_.at(lz-lz_min_);

	hutm.start("RDM");
	for(size_t i=0;i<dim;++i)
	{
		auto basis_b = basis + i*Ne_;
		auto basis_e = basis + (i+1)*Ne_;
		int numThreads = omp_get_max_threads();
		#pragma omp parallel for num_threads(numThreads)
		for(size_t j=0;j<mask_dim_;++j)
			rdm_ij(i, j, lz, dimr, rdm, vector, basis, alphalist, betalist);
		/*{
			int lz_l = 0;
			size_t l_key = 0;
			size_t r_key = 0;
			auto mask_b = l_mask_.cbegin() + j*Nl_;
			auto mask_e = l_mask_.cbegin() + (j+1)*Nl_;
			int l_c = 0;
			int r_c = 0;
			Data_dtype mul = 1.0;

			//hutm.start("count_lz");
			for(auto b_it=mask_b;b_it!=mask_e;b_it++)
			{
				lz_l += basis_b[*b_it]%No_;
			}
			//hutm.stop("count_lz");
			if( lz_l!=lz )
				continue;
			
			//hutm.start("make_key_value");
			for(auto b_it=basis_b;b_it!=basis_e&&mul>1e-30;b_it++)
			{
				if(*b_it==basis_b[*mask_b] && mask_b!=mask_e)// left
				{
					//lz_l += (*b_it)%No_;
					l_key += binomial(*b_it, l_c+1);
					mul *= alphalist[(*b_it)%No_];
					++l_c;
					++mask_b; 
				}
				else// right
				{
					r_key += binomial(*b_it, r_c+1);
					mul *= betalist[(*b_it)%No_];
					++r_c;
				}
			}
			//hutm.stop("make_key_value");
			// add to unordered_map
			auto lit = lhash.find(l_key);
			if(lit==lhash.end())
				continue;
			auto rit = rhash.find(r_key);
			if(rit==rhash.end())
				continue;
			rdm.at(lit->second.ind*dimr+rit->second.ind) = vector[i]*mul*sign_mask_[j];//*lit->second.coef*rit->second.coef;// *phase;
		}*/
	}
	hutm.stop("RDM");
	hutm.start("SVD");
	std::vector<Data_dtype> s(std::min(diml,dimr), 0);
	auto info = svd(rdm.data(), s.data(), diml, dimr);
	hutm.stop("SVD");
	if(info!=0)
		printf("Unconverged! Info=%d\n", info);
	double entropy = 0.0;
	double norm_2 = 0.0;
	for(auto i : s)
	{
		//printf("%lf\n", i);
		if(i>1e-20)
			entropy -= 2.0*i*i*std::log(i);
		norm_2 += i*i;
	}
	/*for(size_t i=0;i<diml;++i)
	{
		for(size_t j=0;j<dimr;++j)
		{
			printf("%.2e ", rdm[i*diml+j]);
		}
		printf("\n");
	}*/
	return std::make_tuple(norm_2, entropy);
}

template<class RandomAccessIterator>
void bipart::
rdm_ij(size_t i, size_t j, int lz, size_t dimr, std::vector<Data_dtype>& rdm, Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist) const
{
	int lz_l = 0;
	size_t l_key = 0;
	size_t r_key = 0;
	auto mask_b = l_mask_.cbegin() + j*Nl_;
	auto mask_e = l_mask_.cbegin() + (j+1)*Nl_;
	auto basis_b = basis + i*Ne_;
	auto basis_e = basis + (i+1)*Ne_;
	int l_c = 0;
	int r_c = 0;
	Data_dtype mul = 1.0;

	//hutm.start("count_lz");
	for(auto b_it=mask_b;b_it!=mask_e;b_it++)
	{
		lz_l += basis_b[*b_it]%No_;
	}
	//hutm.stop("count_lz");
	if( lz_l!=lz )
		return;
	
	const lin_hash& lhash = l_inf_.at(lz-lz_min_);
	const lin_hash& rhash = r_inf_.at(lz-lz_min_);
	//hutm.start("make_key_value");
	for(auto b_it=basis_b;b_it!=basis_e&&mul>1e-30;b_it++)
	{
		if(*b_it==basis_b[*mask_b] && mask_b!=mask_e)// left
		{
			//lz_l += (*b_it)%No_;
			l_key += binomial(*b_it, l_c+1);
			mul *= alphalist[(*b_it)%No_];
			++l_c;
			++mask_b; 
		}
		else// right
		{
			r_key += binomial(*b_it, r_c+1);
			mul *= betalist[(*b_it)%No_];
			++r_c;
		}
	}
	//hutm.stop("make_key_value");
	// add to unordered_map
	auto lit = lhash.find(l_key);
	if(lit==lhash.end())
		return;
	auto rit = rhash.find(r_key);
	if(rit==rhash.end())
		return;
	rdm[lit->second.ind*dimr+rit->second.ind] = vector[i]*mul*sign_mask_[j];//*lit->second.coef*rit->second.coef;// *phase;
}


template<class RandomAccessIterator>
std::tuple<double, double> bipart::
entropy(Data_dtype* vector, RandomAccessIterator basis, Data_dtype* alphalist, Data_dtype* betalist, size_t dim) const
{
	hu::timer hutm;
	hutm.start("end");
	Data_dtype* coef;
	if(Nl_==0)
		coef = betalist;
	else if(Nl_==Ne_)
		coef = alphalist;
	else
		return std::make_tuple(0.0, 0.0);

	int numThreads = omp_get_max_threads();
	Data_dtype result = 0.0;
	#pragma omp parallel for num_threads(numThreads) reduction(+:result)
	for(size_t i=0;i<dim;++i)
	{
		Data_dtype temp = vector[i];
		for(size_t j=0;j<Ne_;++j)
		{
			temp *= coef[basis[i*Ne_+j]%No_];
		}
		result = result + temp*temp;
	}
	double entropy;
	if(result>1e-30)
		entropy = -result*std::log(result);
	else
		entropy = 0.0;
	hutm.stop("end");
	return std::make_tuple(result, entropy);
}

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
		printf("Unconverged! Info=%d\n", info);
	
	mkl_free(superb);
	return info;
}

std::tuple<double, double> 
res_ED_ver1(int N_o, int N_e, std::vector<Data_dtype> vec, std::vector<Data_dtype> & alphalist, std::vector<Data_dtype> & betalist, const struct FQHE& A_)
{
	hu::timer hutm;
	double SA = 0;
	double norm_SA = 0;
	for(auto Nl=0;Nl<=N_e;++Nl)
	{
		hutm.start("bipart");
		auto bi = bipart(N_o, N_e, Nl, -16, A_.basis, A_.dim);
		hutm.stop("bipart");
		hutm.start("entropy");
		for(auto lz=bi.lz_min();lz<=bi.lz_max();lz++)
		{
			auto [norm_2, entropy] = bi.entropy(vec.data(), A_.basis, alphalist.data(), betalist.data(), lz, A_.dim);
			norm_SA += norm_2;
			SA += entropy;
			//printf("\tlz=%d, norm_2 = %lf, entropy = %lf\n", lz, norm_2, entropy);
		}
		hutm.stop("entropy");
		printf("Nl=%d, norm_SA = %lf, entropy = %lf\n", Nl, norm_SA, SA);
	}
	SA = SA/norm_SA+std::log(norm_SA);
	printf("Final EE: S_E/norm+log(norm) = %lf\n", SA);
	return std::make_tuple(SA, 1-norm_SA);
}

std::tuple<double, double, double, double> 
mutual_ED_ver1(int N_o, int N_e, double theta_p, std::vector<Data_dtype> vec, const struct FQHE& A_)
{
	hu::timer hutm;
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

	auto [SA, errorA] = res_ED_ver1(N_o, N_e, vec, alphalist, betalist, A_);

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
	
	auto [SAB, errorAB] = res_ED_ver1(N_o, N_e, vec, alphalist, betalist, A_);

	return std::make_tuple(SA, errorA, SAB, errorAB);
}