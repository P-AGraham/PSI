
void cal2BodyAlist(struct CG2 *cg, double *pseudo, double *A_list2, int num_thread)
{
	double *cgt = NULL;
	long int u = (size_t)cg->u;
	memset(A_list2, 0, u*u*u*u*sizeof(double));

	//printCG2(cg);
	# pragma omp parallel for num_threads(num_thread)
	for(long int m1=0;m1<u;++m1)
		for(long int m2=0;m2<u;++m2)
			for(long int m3=0;m3<u;++m3)
				for(long int m4=0;m4<u;++m4)
					for(long int j=0;j<u;++j)
						A_list2[m1*u*u*u + m2*u*u + m3*u + m4] += cg->cg[u-j-1][(u-1-m1)*u+m2]*cg->cg[u-j-1][(u-1-m4)*u+m3]*pseudo[j];
	//printf("%lf %lf\n",cg->cg[8][(u-1-5)*u+6],cg->cg[8][(u-1-6)*u+5]);
	return;
}

/*

cg table ordering 
			************************************************************************************
			****************************************j=0*****************************************
			************************************************************************************
	3/2 	-5.000000e-01   0.000000e+00    0.000000e+00    0.000000e+00
	1/2 	0.000000e+00    5.000000e-01    0.000000e+00    0.000000e+00
↑	-1/2	0.000000e+00    0.000000e+00    -5.000000e-01   0.000000e+00
m1	-3/2	0.000000e+00    0.000000e+00    0.000000e+00    5.000000e-01  
   m2 ->    -3/2  			-1/2 			1/2 			3/2

usage: <s,m1;s,m2|s,s;j,m1+m2> = cg->cg[u-j-1][(u-1-m1)*u+m2]

*/



void calO00Alist(struct CG2 *cg, double *olist, double *A_list2, int num_thread)
{
	double *cgt = NULL;
	long int u = (size_t)cg->u;
	memset(A_list2, 0, u*u*(2*u-1)*sizeof(double));

	//printCG2(cg);

	# pragma omp parallel for num_threads(num_thread)
	for(long int m1=0;m1<u;++m1)// -s <= m1 <= s
		for(long int m2=0;m2<u;++m2)// -s <= m2 <= s
			for(long int l=0;l<u;++l)// 0 <= l <= 2s
			{
				if( std::abs(olist[l])<1e-15 )
					continue;
				double comm = olist[l] * u*u * cg->cg[l][0] * cg->cg[l][0]/(4*M_PI*std::sqrt(2*l+1));
				// o_l * (2s+1) * <s,s;s,-s|s,s;l,0>^2
				for(long int m=-l;m<l+1;++m)// -l <= m <= l
				{
					if( m1+m<0 || m1+m>=u || m2-m<0 || m2-m>=u )
						continue;
					double dig = comm;
					dig *= ((l+m+m1+m2)&1)?-1:1; // m1-(u-1)/2 + m2-(u-1)/2 = m1+m2-(u-1)
					dig *= cg->cg[l][(u-1-(m1+m))*u+(u-1-m1)] * cg->cg[l][(u-1-(m2-m))*u+(u-1-m2)];
					// (-1)^(l+m+2s+m1+m2) * <s,m1-m;s,-m1|s,s;l,-m> * <s,m2-m;s,-m2|s,s;l,-m>
					//if(l==0)
					//	printf("comm=%lf, cg->cg[%d][%d]=%lf, cg->cg[%d][%d]=%lf\n", 
					//		comm, l, (u-1-(m1+m))*u+(u-1-m1), cg->cg[l][(u-1-(m1+m))*u+(u-1-m1)], l, (u-1-(m2+m))*u+(u-1-m2), cg->cg[l][(u-1-(m2+m))*u+(u-1-m2)]);
					A_list2[m1*u*(2*u-1)+m2*(2*u-1)+(m+u-1)] += dig;
				}
			}
	/*for(long int m1=0;m1<u;++m1)// -s <= m1 <= s
	{
		for(long int m2=0;m2<u;++m2)// -s <= m2 <= s
			printf("%lf\t", A_list2[m1*u*(2*u-1)+m2*2*u+(u-1)]);
		printf("\n");
	}*/
	return;
}


void calOL0Alist(struct CG2 *cg, double *olist, double *A_list2, int L, int num_thread)
{
	double *cgt = NULL;
	long int u = (size_t)cg->u;
	memset(A_list2, 0, u*u*(2*u-1)*sizeof(double));

	struct CG2* cglist = (struct CG2*)malloc( u * sizeof(struct CG2));
	for(int i=0;i<u;++i)
		cglist[i].cg = NULL;
	for(int i=0;i<u;++i)
	{
		printf("%d\n", i);
		int st = 0;
		st = createCG2Bdoy(cglist+i, 2*i+1);
		if(st!=0)
		{
			printf("calOL0Alist::createCG2Bdoy:i=%d\tst=%d\n", i, st);
			goto finish;
		}
		//if(i==0)
		//	printCG2(cglist+i);
		//if(2*i>=L)
		//{
		//	for(long int m=-i;m<i+1;++m)
		//		printf("<l1=%d,m1=%d;l2=%d,m2=%d|L=%d,m=0>=%lf\n", i, -m, i, m, L, cglist[i].cg[L][(m+i)*(2*i+1)+m+i]);
		//}

	}


	# pragma omp parallel for num_threads(num_thread)
	for(long int m1=0;m1<u;++m1)// -s <= m1 <= s
		for(long int m2=0;m2<u;++m2)// -s <= m2 <= s
			for(long int l=0;l<u;++l)// 0 <= l <= 2s
			{
				if( std::abs(olist[l])<1e-15 || 2*l<L )
					continue;
				double comm = olist[l] * u*u * cg->cg[l][0]*cg->cg[l][0]/(4*M_PI);
				// o_l * (2s+1) * <s,s;s,-s|s,s;l,0>^2
				for(long int m=-l;m<l+1;++m)// -l <= m <= l
				{
					if( m1+m<0 || m1+m>=u || m2-m<0 || m2-m>=u )
						continue;
					double dig = comm;
					dig *= ((m1+m2)&1)?-1:1; // m1-(u-1)/2 + m2-(u-1)/2 = m1+m2-(u-1)
					dig *= cglist[l].cg[L][(m+l)*(2*l+1)+m+l]*cg->cg[l][(u-1-(m1+m))*u+(u-1-m1)] * cg->cg[l][(u-1-(m2-m))*u+(u-1-m2)];
					// (-1)^(l+m+2s+m1+m2) * <s,m1-m;s,-m1|s,s;l,-m> * <s,m2-m;s,-m2|s,s;l,-m>
					//if(l==0)
					//	printf("comm=%lf, cg->cg[%d][%d]=%lf, cg->cg[%d][%d]=%lf\n", 
					//		comm, l, (u-1-(m1+m))*u+(u-1-m1), cg->cg[l][(u-1-(m1+m))*u+(u-1-m1)], l, (u-1-(m2+m))*u+(u-1-m2), cg->cg[l][(u-1-(m2+m))*u+(u-1-m2)]);
					A_list2[m1*u*(2*u-1)+m2*(2*u-1)+(m+u-1)] += dig;
				}
			}
	/*for(long int m1=0;m1<u;++m1)// -s <= m1 <= s
	{
		for(long int m2=0;m2<u;++m2)// -s <= m2 <= s
			printf("%lf\t", A_list2[m1*u*(2*u-1)+m2*2*u+(u-1)]);
		printf("\n");
	}*/
	finish:		
	for(int i=0;i<u;++i)
	{
		//printf("%d\n", i);
		destoryCG2(cglist+i);
	}
	free(cglist);
	return;
}


void calnnl0Alist(double *A_list2, int l, int N_o, int num_thread)
{
	double s = 0.5*N_o-0.5;
	memset(A_list2, 0, N_o*N_o*(2*N_o-1)*sizeof(double));
	double comm = N_o*N_o/std::sqrt( (2*l+1)*64*M_PI*M_PI*M_PI );

	auto cglsitss = std::move(cgj1j2(N_o, N_o));

	//# pragma omp parallel for num_threads(num_thread)
	for(int xl1=0;xl1<N_o;++xl1)
	{
		double l1 = (double)xl1;
		for(double L=-l;L<l+0.1;++L)
		{
			if( l1+L<-0.1 || l1+L>2*s+0.1 || l<std::abs(L) || l>std::abs(2*l1+L)+0.1 )
				continue;
			auto cgl1Ll = std::move(cgj1j2j(round(2*l1+1),round(2*l1+2*L+1),round(2*l+1)));
			double temp1 = cglsitss.cgs(l1).cgm(s, -s) * cglsitss.cgs(l1+L).cgm(s, -s) * cgl1Ll.cgm((double)0.0, (double)0.0);
			for(double m=-std::min(l1,l1+L); m<std::min(l1,l1+L)+0.1; ++m)
			{
				double temp2 = temp1 * cgl1Ll.cgm(-m, m);
				for(double m1=-s;m1<s+0.1;++m1)
				{
					for(double m2=-s;m2<s+0.1;++m2)
					{
						if( std::abs(m1+m)>s+0.1 || std::abs(m2-m)>s+0.1 )
							continue;
						double dig = temp2 * cglsitss.cgs(l1).cgm(m1+m, -m1) * cglsitss.cgs(l1+L).cgm(m2-m, -m2) * comm;
						dig *= ((int)round(std::abs(2*s+2*l1+L+m1+m2))&1)?-1:1;
						A_list2[(int)round(m1+s)*N_o*(2*N_o-1)+(int)round(m2+s)*(2*N_o-1)+(int)round(m+N_o-1)] += dig;
						//A_list2[m1*u*(2*u-1)+m2*(2*u-1)+(m+u-1)] += dig;
					}
				}
			}
		}
	}
}

// sum_m n_{l,-m} n_{l,m}
void cal_nlm_nlm_Alist(double *A_list2, double *olist, int l, int N_o, int num_thread)
{
	double s = 0.5*N_o-0.5;
	memset(A_list2, 0, N_o*N_o*(2*N_o-1)*sizeof(double));
	double comm = N_o*N_o/( 4*M_PI*(2*l+1) );

	auto cg = std::move(cgj1j2j(N_o, N_o, round(2*l+1)));

	comm *= cg.cgm(s,-s) * cg.cgm(s,-s);

	for(double m=-l;m<l+0.1;++m)
	{
		for(double m1=-s;m1<s+0.1;++m1)
		{
			for(double m2=-s;m2<s+0.1;++m2)
			{
				if( std::abs(m1+m)>s+0.1 || std::abs(m2-m)>s+0.1 )
					continue;
				double dig = cg.cgm(m1+m,-m1) * cg.cgm(m2-m,-m2) * comm * olist[(int)round(m+l)];
				dig *= ((int)round(std::abs(2*s+m1+m2))&1)?-1:1;
				A_list2[(int)round(m1+s)*N_o*(2*N_o-1)+(int)round(m2+s)*(2*N_o-1)+(int)round(m+N_o-1)] += dig;
			}
		}		
	}
}


// n_{l_1,-m} n_{l_2,m}
void cal_nl1m_nl2m_Alist(double *A_list2, int l_1, int l_2, int m, int N_o, int num_thread)
{
	//printf("l1=%d,l2=%d\n", l_1, l_2);
	double s = 0.5*N_o-0.5;
	memset(A_list2, 0, N_o*N_o*(2*N_o-1)*sizeof(double));
	double comm = N_o*N_o/( 4*M_PI*std::sqrt((2*l_1+1)*(2*l_2+1)) );
	comm *= ((N_o-1+l_1+l_2)&1)?-1:1;

	auto cg1 = std::move(cgj1j2j(N_o, N_o, 2*l_1+1));
	auto cg2 = std::move(cgj1j2j(N_o, N_o, 2*l_2+1));
	//std::cout<< cg1 << std::endl;
	//std::cout<< cg2 << std::endl;

	comm *= cg1.cgm(s,-s) * cg2.cgm(s,-s);

	for(double m1=-s;m1<s+0.1;++m1)
	{
		//printf("m1=%lf, comm=%lf\n", m1, comm);
		for(double m2=-s;m2<s+0.1;++m2)
		{
			if( std::abs(m1+m)>s+0.1 || std::abs(m2-m)>s+0.1 )
				continue;
			double dig = cg1.cgm(m1+m,-m1) * cg2.cgm(m2-m,-m2) * comm;
			dig *= ((int)round(std::abs(m1+m2))&1)?-1:1;
			//printf("\tm1=%lf,m2=%lf, cg1.cgm(%lf, %lf)=%lf, cg1.cgm(%lf, %lf)=%lf, ,resu=%lf\n", m1, m2, m1+m, -m1, cg1.cgm(m1+m,-m1), m2-m, -m2, cg2.cgm(m2-m,-m2), dig);
			A_list2[(int)round(m1+s)*N_o*(2*N_o-1)+(int)round(m2+s)*(2*N_o-1)+(int)round(m+N_o-1)] += dig;
		}		
	}
}


// pairing operator \Delta
void cal_delta_Alist(double *A_list2, int l, int m, int N_o, int num_thread)
{
	//printf("l1=%d,l2=%d\n", l_1, l_2);
	double s = 0.5*N_o-0.5;
	memset(A_list2, 0, N_o*N_o*sizeof(double));
	//double comm = N_o*N_o/( 4*M_PI*std::sqrt((2*l_1+1)*(2*l_2+1)) );
	//comm *= ((N_o-1+l_1+l_2)&1)?-1:1;

	auto cg1 = std::move(cgj1j2j(N_o, N_o, 2*l+1));

	for(double m1=-s;m1<s+0.1;++m1)
	{
		//printf("m1=%lf, comm=%lf\n", m1, comm);
		if( std::abs(m-m1)>s+0.1 || std::abs(m1)>s+0.1 )
			continue;
		//double dig = cg1.cgm(m1+m,-m1) * cg2.cgm(m2-m,-m2) * comm;
		//dig *= ((int)round(std::abs(m1+m2))&1)?-1:1;
		//printf("\tm1=%lf,m2=%lf, cg1.cgm(%lf, %lf)=%lf, cg1.cgm(%lf, %lf)=%lf, ,resu=%lf\n", m1, m2, m1+m, -m1, cg1.cgm(m1+m,-m1), m2-m, -m2, cg2.cgm(m2-m,-m2), dig);
		A_list2[(int)round(m1+s)*N_o+(int)round(m-m1+s)] += cg1.cgm(m1,m-m1);
	}
}


// pairing operator \Delta*\Delta
void cal_deltadelta_lalb_Alist(double *A_list2, int l, int m, int la, int lb, int N_o, int num_thread)
{
	double s = 0.5*N_o-0.5;
	memset(A_list2, 0, N_o*N_o*N_o*N_o*sizeof(double));
	/*// my suggestion
	double comm = (std::abs(la-lb)&1)?-1:1;
	comm /= std::sqrt(2*l+1);*/
	// Prof. He's suggestion
	double comm = (std::abs(la-lb+m)&1)?-1:1;
	comm /= std::sqrt(2*l+1);

	auto cglalbl = std::move(cgj1j2j( round(2*la+1), round(2*lb+1), round(2*l+1) ));
	auto cgssla = std::move(cgj1j2j(N_o, N_o, round(2*la+1) ));
	auto cgsslb = std::move(cgj1j2j(N_o, N_o, round(2*lb+1) ));

	//# pragma omp parallel for num_threads(num_thread)
	for(double j1=-s;j1<s+0.1;++j1)
	{
		for(double j2=-s;j2<s+0.1;++j2)
		{
			for(double j3=-s;j3<s+0.1;++j3)
			{
				for(double j4=-s;j4<s+0.1;++j4)
				{
					if( std::abs(j1+j2)>la+0.1 || std::abs(j3+j4)>lb+0.1 )
						continue;
					if( std::abs(j1+j2-j3-j4-m)>1e-14 )
						continue;
					A_list2[(int)round(j1+s)*N_o*N_o*N_o+(int)round(j2+s)*N_o*N_o+(int)round(j3+s)*N_o+(int)round(j4+s)] = 
					comm * cglalbl.cgm(j1+j2,-j3-j4) * cgssla.cgm(j1, j2) * cgsslb.cgm(j4, j3);
					// Prof. He's suggestion
					A_list2[(int)round(j1+s)*N_o*N_o*N_o+(int)round(j2+s)*N_o*N_o+(int)round(j3+s)*N_o+(int)round(j4+s)] *= 
					((int)round(j3+j4)&1)?-1:1;
				}
			}
		}		
	}

	/*#ifdef __PRINTCG
	printf("l=%d, la=%d, lb=%d\n", l ,la ,lb);
	std::cout << cglalbl << std::endl;
	std::cout << cgssla  << std::endl;
	std::cout << cgsslb  << std::endl;
	#endif
	for(double j1=-s;j1<s+0.1;++j1)
	{
		for(double j2=-s;j2<s+0.1;++j2)
		{
			for(double j3=-s;j3<s+0.1;++j3)
			{
				for(double j4=-s;j4<s+0.1;++j4)
				{
					if( std::abs(j1+j2)>la+0.1 || std::abs(j3+j4)>lb+0.1 )
						continue;
					if( std::abs(j1+j2-j3-j4-m)>1e-14 )
						continue;
					printf("%d, %d, %d, %d:\t%lf*%lf*%lf*%lf = %lf\n", (int)round(j1+s), (int)round(j2+s), (int)round(j3+s), (int)round(j4+s), 
							comm, cglalbl.cgm(j1+j2,-j3-j4), cgssla.cgm(j1, j2), cgsslb.cgm(j4, j3),
							A_list2[(int)round(j1+s)*N_o*N_o*N_o+(int)round(j2+s)*N_o*N_o+(int)round(j3+s)*N_o+(int)round(j4+s)]);
				}
			}
		}		
	}*/
}

// nlm operator
void cal_nlm_Alist(double *A_list2, int l, int m, int N_o, int num_thread)
{
	double s = 0.5*N_o-0.5;
	memset(A_list2, 0, N_o*N_o*sizeof(double));

	double comm = (std::abs(N_o-1+l+m)&1)?-1:1; // 3s+l+m+j1 = (2s+l+m)+(s+j1)
	comm *= std::sqrt( N_o*N_o/(4*M_PI*(2*l+1)) );

	auto cgssl = std::move(cgj1j2j( N_o, N_o, round(2*l+1) ));

	//# pragma omp parallel for num_threads(num_thread)
	for(double j1=-s;j1<s+0.1;++j1)
	{
		for(double j2=-s;j2<s+0.1;++j2)
		{
			if( std::abs(j1-j2-m)>1e-14 )
				continue;
			A_list2[(int)round(j1+s)*N_o+(int)round(j2+s)] = comm * cgssl.cgm(s,-s) * cgssl.cgm(j2,-j1) * (((int)round(j1+s)&1)?-1:1);
		}		
	}

	/*#ifdef __PRINTCG
	printf("l=%d, m=%d\n", l, m);
	std::cout << cgssl << std::endl;
	#endif

	for(double j1=-s;j1<s+0.1;++j1)
	{
		for(double j2=-s;j2<s+0.1;++j2)
		{
			if( std::abs(j1-j2-m)>1e-14 )
				continue;
			printf("%d, %d:\t%lf*%lf*%lf = %lf\n", (int)round(j1+s), (int)round(j2+s), 
					comm, cgssl.cgm(s,-s), cgssl.cgm(j2,-j1), A_list2[(int)round(j1+s)*N_o+(int)round(j2+s)]);
		}
	}*/
}