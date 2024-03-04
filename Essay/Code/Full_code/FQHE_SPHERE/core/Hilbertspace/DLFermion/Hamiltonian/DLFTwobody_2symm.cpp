

void basis_symm::element_count_thread(unsigned int row, unsigned int* indptr)
{
	unsigned int bra_ind, sub_ind, bra_num, ind, n_data = 0;
	int j1, j2, j3, j4, J12, J34, phase, phase1, phase2;
	size_t component_num;
	const unsigned char *basis_t, *basis_u, *basis_d;

	unsigned char* bra_i = new unsigned char [N_phi];
	unsigned char* bra_j = new unsigned char [N_phi]; //N_e
	unsigned char* diff_u = new unsigned char [N_phi];
	unsigned char* diff_d = new unsigned char [N_phi];

	unordered_map<size_t, size_t>::iterator hash_it;
	unordered_map<size_t, size_t>::iterator hash_ed = k2_hash.end();
	basis *temp = k1_basis + k2_basis[row];
	//cout << row <<"\t"<< indptr[(size_t)row] << "\t";temp->print(N_u,N_d);cout << endl;
	//----------------------------------------------------------------------------------------------
	// the first component
	//-----------------------------------------------------------------------------------------------
	size_t k;
	if(N_u>0)
	{
		k = 0;
		basis_t = temp->up();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_u))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_u[k] = j;
				++k;
			}
		}
		basis_t = temp->down();
		component_num = dict_num(basis_t, N_d);
		basis_t = temp->up();

		for(size_t u=0;u<N_u-1;++u)
		{
			j3 = basis_t[u];
			for(size_t v=u+1;v<N_u;++v)
			{
				j4 = basis_t[v];
				J34 = j3 + j4;
				if(J34>=N_phi)
					J34 -= N_phi;
				for(size_t s=0;s<N_phi-N_u-1;++s)
				{
					j2 = diff_u[s];
					for(size_t t=s+1;t<N_phi-N_u;++t)
					{
						j1 = diff_u[t];
						J12 = j1 + j2;
						if(J12>=N_phi)
							J12 -= N_phi;
						if(J12!=J34)
							continue;
						//cout<<"2";getchar();
						bra_num = hopping(basis_t, bra_i, bra_j, N_phi, N_u, j1, j2, j3, j4, &phase);
						bra_num = bra_num * stride + component_num;
						ind = find(bra_num);
						bra_ind = k1_basis[ind].first;
						hash_it = k2_hash.find(bra_ind);
						//printf("\t1st\t%u:%d %d %d %d %u %u %u\n",row,j1,j2,j3,j4,ind,bra_ind,hash_it->second);
						//cout<<"3";getchar();
						if(hash_it!=hash_ed)
							++n_data;
					}
				}
			}
		}
		++n_data;
	}

	//----------------------------------------------------------------------------------------------
	// the second component
	//-----------------------------------------------------------------------------------------------

	if(N_d>0)
	{
		k = 0;
		basis_t = temp->down();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_d))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_d[k] = j;
				++k;
			}
		}
		basis_t = temp->up();
		component_num = dict_num(basis_t, N_u);
		basis_t = temp->down();

		for(size_t u=0;u<N_d-1;++u)
		{
			j3 = basis_t[u];
			for(size_t v=u+1;v<N_d;++v)
			{
				j4 = basis_t[v];
				J34 = j3 + j4;
				if(J34>=N_phi)
					J34 -= N_phi;
				for(size_t s=0;s<N_phi-N_d-1;++s)
				{
					j2 = diff_d[s];
					for(size_t t=s+1;t<N_phi-N_d;++t)
					{
						j1 = diff_d[t];
						J12 = j1 + j2;
						if(J12>=N_phi)
							J12 -= N_phi;
						if(J12!=J34)
							continue;
						bra_num = hopping(basis_t, bra_i, bra_j, N_phi, N_d, j1, j2, j3, j4, &phase);
						bra_num = bra_num + component_num * stride;
						ind = find(bra_num);
						bra_ind = k1_basis[ind].first;
						hash_it = k2_hash.find(bra_ind);
						//printf("\t2nd\t%u:%d %d %d %d %u %u %u\n",row,j1,j2,j3,j4,ind,bra_ind,hash_it->second);
						//cout<<"4";getchar();
						if(hash_it!=hash_ed)
							++n_data;
					}
				}
			}
		}
		++n_data;
	}


	//----------------------------------------------------------------------------------------------
	// the interation of two components
	//-----------------------------------------------------------------------------------------------

	if(N_u>0 && N_d>0)
	{
		k = 0;
		basis_t = temp->up();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_u))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_u[k] = j;
				++k;
			}
		}
		k = 0;
		basis_t = temp->down();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_d))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_d[k] = j;
				++k;
			}
		}
		basis_u = temp->up();
		basis_d = temp->down();

		for(size_t u=0;u<N_d;++u)
		{
			j3 = basis_d[u];
			for(size_t v=0;v<N_u;++v)
			{
				j4 = basis_u[v];
				J34 = j3 + j4;
				if(J34>=N_phi)
					J34 -= N_phi;
				for(size_t s=0;s<N_phi-N_d;++s)
				{
					j2 = diff_d[s];
					for(size_t t=0;t<N_phi-N_u;++t)
					{
						j1 = diff_u[t];
						J12 = j1 + j2;
						if(J12>=N_phi)
							J12 -= N_phi;
						if(J12!=J34)
							continue;
						bra_num = hoppingOneBody(basis_u, bra_i, bra_j, N_phi, N_u, j1, j4, &phase1);
						bra_num = bra_num * stride + hoppingOneBody(basis_d, bra_i, bra_j, N_phi, N_d, j2, j3, &phase2);
						ind = find(bra_num);
						bra_ind = k1_basis[ind].first;
						hash_it = k2_hash.find(bra_ind);
						//printf("\t12\t%u:%d %d %d %d %u %u %u\n",row,j1,j2,j3,j4,ind,bra_ind,hash_it->second);
						//cout<<"5";getchar();
						if(hash_it!=hash_ed)
							++n_data;
					}
				}
			}
		}
		++n_data;
	}
	indptr[row+1] = n_data;
	delete [] bra_i;//free(bra_i);
	delete [] bra_j;//free(bra_j);
	delete [] diff_u;//free(diff_u);
	delete [] diff_d;//free(diff_d);
	return;
}

size_t basis_symm::element_count(unsigned int* indptr, int num_thread)
{
	unsigned int i;
	indptr[0] = 0;

	# pragma omp parallel for num_threads(num_thread)
	for(i=0;i<k2_dim;++i)
		element_count_thread(i, indptr);

	for(i=1;i<=k2_dim;++i)
		indptr[i] += indptr[i-1];
	return indptr[k2_dim];
}

size_t basis_symm::element_generate_thread(unsigned int row, \
			complex<double>* A1_list, complex<double>* A2_list, complex<double>* A12_list,\
			unsigned int* indptr, unsigned int* indices, complex<double>* data)
{
	unsigned int bra_ind, sub_ind, bra_num, ind, n_data = 0;
	int j1, j2, j3, j4, j13, j14, j23, j24, J12, J34, phase, phase1, phase2;
	size_t component_num;
	const unsigned char *basis_t, *basis_u, *basis_d;
	unsigned int* indices_thread = indices + (size_t)indptr[row];
	complex<double>* data_thread = data + (size_t)indptr[row];

	unsigned char* bra_i = new unsigned char [N_phi];
	unsigned char* bra_j = new unsigned char [N_phi]; //N_e
	unsigned char* diff_u = new unsigned char [N_phi];
	unsigned char* diff_d = new unsigned char [N_phi];

	unordered_map<size_t, size_t>::iterator hash_it;
	unordered_map<size_t, size_t>::iterator hash_ed = k2_hash.end();
	complex<double> dig = 0.0;

	basis *temp = k1_basis + k2_basis[row];
	//----------------------------------------------------------------------------------------------
	// the first component
	//-----------------------------------------------------------------------------------------------
	size_t k;
	if(N_u>0)
	{
		k = 0;
		basis_t = temp->up();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_u))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_u[k] = j;
				++k;
			}
		}
		basis_t = temp->down();
		component_num = dict_num(basis_t, N_d);
		basis_t = temp->up();

		for(size_t u=0;u<N_u-1;++u)
		{
			j3 = basis_t[u];
			for(size_t v=u+1;v<N_u;++v)
			{
				j4 = basis_t[v];
				J34 = j3 + j4;
				if(J34>=N_phi)
					J34 -= N_phi;
				for(size_t s=0;s<N_phi-N_u-1;++s)
				{
					j2 = diff_u[s];
					for(size_t t=s+1;t<N_phi-N_u;++t)
					{
						j1 = diff_u[t];
						J12 = j1 + j2;
						if(J12>=N_phi)
							J12 -= N_phi;
						if(J12!=J34)
							continue;
						//cout<<"2";getchar();
						bra_num = hopping(basis_t, bra_i, bra_j, N_phi, N_u, j1, j2, j3, j4, &phase);
						bra_num = bra_num * stride + component_num;
						ind = find(bra_num);
						bra_ind = k1_basis[ind].first;
						hash_it = k2_hash.find(bra_ind);
						//printf("%u:%d %d %d %d %u %u %u\n",row,j1,j2,j3,j4,ind,bra_ind,hash_it->second);
						//cout<<"3";getchar();
						if(hash_it!=hash_ed)
						{
							indices_thread[n_data] = hash_it->second;
							j13 = N_phi-1+j1-j3;
							j14 = j13+j3-j4;
							j23 = j13-j1+j2;
							j24 = j14-j1+j2;
							dig = 2.0 * (A1_list[j13*(2*N_phi-1)+j14]-A1_list[j23*(2*N_phi-1)+j24]);
							dig *= sqrt((double)(temp->period)/(k1_basis[bra_ind].period));
							dig *= exp(-2.0*ii*M_PI*k2*k1_basis[ind].sn/pqgcd);
							data_thread[n_data] = dig * phase;
							if(k1_basis[ind].sn!=0)
								data_thread[n_data] *= k1_basis[ind].phase;
							++n_data;
						}
					}
				}
			}
		}
		dig = 0.0;
		for(size_t s1=0;s1<N_u-1;++s1)
		{
			for(size_t s2=s1+1;s2<N_u;++s2)
			{
				j3 = basis_t[s2];
				j4 = basis_t[s1];
				j1 = j4;
				j2 = j3;
				j13 = N_phi-1+j1-j3;
				j14 = j13+j3-j4;
				j23 = j13-j1+j2;
				j24 = j14-j1+j2;
				dig += 2.0 * (A1_list[j13*(2*N_phi-1)+j14]-A1_list[j23*(2*N_phi-1)+j24]);
			}
		}
		indices_thread[n_data] = row;
		data_thread[n_data] = dig;
		++n_data;
	}

	//----------------------------------------------------------------------------------------------
	// the second component
	//-----------------------------------------------------------------------------------------------
	if(N_d>0)
	{
		k = 0;
		basis_t = temp->down();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_d))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_d[k] = j;
				++k;
			}
		}
		basis_t = temp->up();
		component_num = dict_num(basis_t, N_u);
		basis_t = temp->down();

		for(size_t u=0;u<N_d-1;++u)
		{
			j3 = basis_t[u];
			for(size_t v=u+1;v<N_d;++v)
			{
				j4 = basis_t[v];
				J34 = j3 + j4;
				if(J34>=N_phi)
					J34 -= N_phi;
				for(size_t s=0;s<N_phi-N_d-1;++s)
				{
					j2 = diff_d[s];
					for(size_t t=s+1;t<N_phi-N_d;++t)
					{
						j1 = diff_d[t];
						J12 = j1 + j2;
						if(J12>=N_phi)
							J12 -= N_phi;
						if(J12!=J34)
							continue;
						bra_num = hopping(basis_t, bra_i, bra_j, N_phi, N_d, j1, j2, j3, j4, &phase);
						bra_num = bra_num + component_num * stride;
						ind = find(bra_num);
						bra_ind = k1_basis[ind].first;
						hash_it = k2_hash.find(bra_ind);
						//printf("%u:%d %d %d %d %u %u %u\n",row,j1,j2,j3,j4,ind,bra_ind,hash_it->second);
						//cout<<"4";getchar();
						if(hash_it!=hash_ed)
						{
							indices_thread[n_data] = hash_it->second;
							j13 = N_phi-1+j1-j3;
							j14 = j13+j3-j4;
							j23 = j13-j1+j2;
							j24 = j14-j1+j2;
							dig = 2.0 * (A2_list[j13*(2*N_phi-1)+j14]-A2_list[j23*(2*N_phi-1)+j24]);
							dig *= sqrt((double)(temp->period)/(k1_basis[bra_ind].period));
							dig *= exp(-2.0*ii*M_PI*k2*k1_basis[ind].sn/pqgcd);
							data_thread[n_data] = dig * phase;
							if(k1_basis[ind].sn!=0)
								data_thread[n_data] *= k1_basis[ind].phase;
							++n_data;
						}
					}
				}
			}
		}
		dig = 0.0;
		for(size_t s1=0;s1<N_d-1;++s1)
		{
			for(size_t s2=s1+1;s2<N_d;++s2)
			{
				j3 = basis_t[s2];
				j4 = basis_t[s1];
				j1 = j4;
				j2 = j3;
				j13 = N_phi-1+j1-j3;
				j14 = j13+j3-j4;
				j23 = j13-j1+j2;
				j24 = j14-j1+j2;
				dig += 2.0 * (A2_list[j13*(2*N_phi-1)+j14]-A2_list[j23*(2*N_phi-1)+j24]);
			}
		}
		indices_thread[n_data] = row;
		data_thread[n_data] = dig;
		++n_data;
	}

	//----------------------------------------------------------------------------------------------
	// the interation of two components
	//-----------------------------------------------------------------------------------------------

	if(N_u>0 && N_d>0)
	{
		k = 0;
		basis_t = temp->up();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_u))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_u[k] = j;
				++k;
			}
		}
		k = 0;
		basis_t = temp->down();
		for(size_t j=0,ll=0;j<N_phi;++j)
		{
			if((*basis_t == j)&&(ll<N_d))
			{
				++basis_t;
				++ll;
			}
			else
			{
				diff_d[k] = j;
				++k;
			}
		}
		basis_u = temp->up();
		basis_d = temp->down();

		for(size_t u=0;u<N_d;++u)
		{
			j3 = basis_d[u];
			for(size_t v=0;v<N_u;++v)
			{
				j4 = basis_u[v];
				J34 = j3 + j4;
				if(J34>=N_phi)
					J34 -= N_phi;
				for(size_t s=0;s<N_phi-N_d;++s)
				{
					j2 = diff_d[s];
					for(size_t t=0;t<N_phi-N_u;++t)
					{
						j1 = diff_u[t];
						J12 = j1 + j2;
						if(J12>=N_phi)
							J12 -= N_phi;
						if(J12!=J34)
							continue;
						bra_num = hoppingOneBody(basis_u, bra_i, bra_j, N_phi, N_u, j1, j4, &phase1);
						bra_num = bra_num * stride + hoppingOneBody(basis_d, bra_i, bra_j, N_phi, N_d, j2, j3, &phase2);
						ind = find(bra_num);
						bra_ind = k1_basis[ind].first;
						hash_it = k2_hash.find(bra_ind);
						//printf("%u:%d %d %d %d %u %u %u\n",row,j1,j2,j3,j4,ind,bra_ind,hash_it->second);
						//cout<<"5";getchar();
						if(hash_it!=hash_ed)
						{
							indices_thread[n_data] = hash_it->second;
							j13 = N_phi-1+j1-j3;
							j14 = j13+j3-j4;
							//printf("%u:%d %d %d %d ",row,j1,j2,j3,j4);
							dig = A12_list[j13*(2*N_phi-1)+j14];//cout << dig << endl;
							dig *= sqrt(((double)(temp->period))/(k1_basis[bra_ind].period));
							dig *= exp(-2.0*ii*M_PI*k2*k1_basis[ind].sn/pqgcd);
							data_thread[n_data] = dig * phase1 * phase2;
							if(k1_basis[ind].sn!=0)
								data_thread[n_data] *= k1_basis[ind].phase;
							++n_data;
						}
					}
				}
			}
		}
		dig = 0.0;
		for(size_t s1=0;s1<N_u;++s1)
		{
			for(size_t s2=0;s2<N_d;++s2)
			{
				j3 = basis_d[s2];
				j4 = basis_u[s1];
				j1 = j4;
				j13 = N_phi-1+j1-j3;
				j14 = j13+j3-j4;
				dig += A12_list[j13*(2*N_phi-1)+j14];
			}
		}
		indices_thread[n_data] = row;
		data_thread[n_data] = dig;
		++n_data;
	}
	//printf("\n");
	delete [] bra_i;//free(bra_i);
	delete [] bra_j;//free(bra_j);
	delete [] diff_u;//free(diff_u);
	delete [] diff_d;//free(diff_d);
	return n_data;
}

size_t basis_symm::element_generate(complex<double>* A1_list, \
			complex<double>* A2_list, complex<double>* A12_list,  unsigned int* indptr, \
			unsigned int* indices, complex<double>* data, int num_thread)
{
	unsigned int i;
	indptr[0] = 0;

	# pragma omp parallel for num_threads(num_thread)
	for(i=0;i<k2_dim;++i)
	{
		if(element_generate_thread(i, A1_list, A2_list, A12_list, \
						indptr, indices, data) != indptr[i+1] - indptr[i])
		{
			cout << "row = " << i << endl;exit(10);
		}
	}

	return indptr[k2_dim];
}
