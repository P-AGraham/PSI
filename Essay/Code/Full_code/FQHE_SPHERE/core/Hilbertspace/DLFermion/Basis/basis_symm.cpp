

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<class basis_symm>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
/*class basis_symm
{
private:
	basis *k1_basis;
	basis *ket;
	size_t N_u;
	size_t N_d;
	size_t N_phi;
	size_t k1_dim;
	size_t k1;
	size_t stride;
	size_t p;
	size_t q;
	size_t pqgcd;
	unordered_map<size_t, size_t> k1_hash;
public:
	basis_symm(){k1_basis=NULL;}
	void create(basis_array &u, basis_array &d, size_t k);
	~basis_symm();

	friend std::ostream & operator<<(std::ostream & os, basis_symm & b);

	void operator++(){++ket;}
	void operator--(){--ket;}
	void set(size_t i){ket=k1_basis+i;}
	size_t find(size_t key){return k1_hash[key];}
	size_t key(){ket->base_num=dict_num(ket->up(),N_u)*stride+dict_num(ket->down(),N_d);return ket->base_num;}
	size_t find(){return k1_hash[this->key()];} // cal the key of *ket and return the value k1_hash[key];
	void translation();
	
	void hash();
};*/

void basis_symm::create_k2(size_t k)
{
	if(k2==k)
		return;
	k2 = k;
	k2_hash.clear();
	if(k2_basis!= NULL)
		delete [] k2_basis;
	size_t j=0, q1;
	for(size_t i=0;i<k1_dim;++i)
	{
		this->set(i);
		//ket->print(N_u, N_d);cout << endl;
		if(ket->period==0)
			continue;
		if(ket->phase==1)//exp(iphi)=1 phi=0 for bilayer Halperin phi = 0
		{
			if((ket->period*k2) % (pqgcd))
				continue;
			k2_hash[i] = j;//cout << j << endl;
			++j;
		}
		else//exp(iphi)=-1 phi=pi
		{
			q1 = (pqgcd)/ket->period;
			if((k2+q1/2)%q1)
				continue;
			k2_hash[i] = j;//cout << j << endl;
			++j;
		}
	}
	k2_dim = j;
	k2_basis = new size_t [k2_dim];
	unordered_map<size_t, size_t>::iterator hash_it;
	unordered_map<size_t, size_t>::iterator hash_ed = k2_hash.end();
	for(size_t i=0;i<k1_dim;++i)
	{
		//this->set(i);
		hash_it = k2_hash.find(i);
		if(hash_it != hash_ed)
		{
			//cout << hash_it->first << " -> " << hash_it->second;
			//cout << " -> " << i << endl;
			k2_basis[hash_it->second] = i;
		}
	}
	this->set(0);
}

void basis_symm::translational_symm()
{
	if(translational_symm_Q!=0)
		return;
	size_t temp, temp1, temp2, j, ind;
	char phase1, phase2;
	for(size_t i=0;i<k1_dim;++i)
	{
		this->set(i);
		if(ket->sn!=0) //!=0
			continue;
		for(j=0;j<pqgcd;++j)
		{
			this->set(i);
			temp1 = base_translation(ket->up(), q*(j+1), N_phi, N_u, &phase1);
			temp2 = base_translation(ket->down(), q*(j+1), N_phi, N_d, &phase2);
			temp = temp1*stride + temp2;
			if(temp == ket->base_num)
				break;
			ind = find(temp);
			this->set(ind);
			ket->first = i;
			ket->sn = j + 1;
			ket->phase = phase1 * phase2;
		}
		this->set(i);
		ket->first = i;
		ket->period = j + 1;
		ket->phase = phase1 * phase2;//(((j+1)*p*(N_e - 1))&1==1)?-1:1;//(((j+1)*p*(N_e - 1))&1==1)?-1:1;
	}
	translational_symm_Q = 1;
	this->set(0);
}


void basis_symm::hash()
{
	for(size_t i=0;i<k1_dim;++i)
	{
		this->set(i);
		this->key();
		k1_hash[ket->base_num] = i;
	}
	this->set(0);
}

void basis_symm::create_k1(basis_array &u, basis_array &d, size_t k)
{
	if(k1_basis!=NULL)
		delete [] k1_basis;
	k1_basis = NULL;
	if(k2_basis!=NULL)
		delete [] k2_basis;
	k2_basis = NULL;
	k1_hash.clear();
	k2_hash.clear();
	translational_symm_Q = 0;
	N_u = u.get_N_e();
	N_d = d.get_N_e();
	N_phi = u.get_N_phi();
	k1 = k;
	k2 = 0xFFFF;
	k1_dim = 0;
	k2_dim = 0;

	size_t dim_u, dim_d;
	unsigned char m, n;
	dim_u = u.get_dim();
	dim_d = d.get_dim();

	for(size_t i=0;i<dim_u;++i)
	{
		m = u.get_k(i);
		for(size_t j=0;j<dim_d;++j)
		{
			n = d.get_k(j);
			if((m+n==k1)||(m+n==k1+N_phi))
				++k1_dim;
		}
	}
	
	size_t count = 0;
	k1_basis = new basis[k1_dim];
	for(size_t i=0;i<dim_u;++i)
	{
		m = u.get_k(i);
		for(size_t j=0;j<dim_d;++j)
		{
			n = d.get_k(j);
			if((m+n==k1)||(m+n==k1+N_phi))
			{
				k1_basis[count].init(u[i], d[j], count);
				++count;
			}
		}
	}
	ket = k1_basis;
	stride = binomial(N_phi, N_d);
	pqgcd = gcd(N_phi, N_u+N_d);
	p = (N_u+N_d)/pqgcd;
	q = N_phi/pqgcd;
}

void basis_symm::reset()
{
	if(k1_basis!=NULL)
		delete [] k1_basis;
	k1_basis = NULL;
	if(k2_basis!=NULL)
		delete [] k2_basis;
	k2_basis = NULL;
	k1_hash.clear();
	k2_hash.clear();
	translational_symm_Q = 0;
	k1 = 0xFFFF;
	k2 = 0xFFFF;
	k1_dim = 0;
	k2_dim = 0;
}


basis_symm::~basis_symm()
{
	//cout << "destroy basis_symm at " << this << endl;
	//cout << "aa" << k1_basis << endl;
	if(k1_basis!=NULL)
		delete [] k1_basis;
	if(k2_basis!=NULL)
		delete [] k2_basis;
	//k1_basis = NULL;
	k1_hash.clear();
	k2_hash.clear();
}

void basis::print(size_t N_up, size_t N_down)
{
	cout << index << "\tu: ";
	for(size_t i=0;i<N_up;++i)
		cout << (int)u[i] << " ";
	cout << "\td: ";
	for(size_t i=0;i<N_down;++i)
		cout << (int)d[i] << " ";
	cout << "\t";
	cout << "index=" << index << " ";
	cout << "first=" << first << " ";
	cout << "period=" << (int)period << " ";
	cout << "sn=" << (int)sn << " ";
	cout << "phase=" << (int)phase << "\t ";
	cout << "hash " << base_num;
}

std::ostream & operator<<(std::ostream & os, basis_symm & b)
{
	unordered_map<size_t, size_t>::iterator hash_it;
	unordered_map<size_t, size_t>::iterator hash_ed = b.k2_hash.end();
	for(size_t i=0;i<b.k1_dim;++i)
	{
		b.set(i);
		hash_it = b.k2_hash.find(i);
		b.k1_basis[i].print(b.N_u, b.N_d);
		cout << " value:" << b.find(b.k1_basis[i].base_num);
		if((hash_it!=hash_ed))
			cout << "\t" << hash_it->second;
		else
			cout << "\t" << -1;
		cout << endl;
	}
	cout << "dim1 = " << b.k1_dim << "\tdim2 = " << b.k2_dim << endl;
	b.set(0);
	return os;
}

void basis_symm::copy_basis(unsigned char *out, size_t i)
{
	for(size_t j=0;j<N_u;++j)
		out[j] = k1_basis[i].up()[j];
	for(size_t j=0;j<N_d;++j)
		out[j+N_u] = k1_basis[i].down()[j];
}