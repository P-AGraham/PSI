

/*
	hu::timer hutm;	

	if(argc < 4) 
	{ 
	   //reminds us to give an input file if we forget
	   printf("Usage: %s 1001001001000 3 1\n",argv[0]); 
	   return 0; 
	}
	std::string root_config(argv[1]);
	int q = std::stoi(argv[2]);
	int p = std::stoi(argv[3]);
	int dk = 5;
	bool count0Q = false;
	if(argc>4)
		dk = std::stoi(argv[4]);

	size_t N_o = root_config.length();
	size_t N_e = 0;//std::count( root_config.begin(), root_config.end(), '1' );
	int K = 0;
	int Kq = 0;
	Data_basis *temp =  (Data_basis *)malloc((size_t)A.N_e * sizeof(Data_basis));
	Data_basis *temp1 = temp;
	for(auto i=0;i<N_o;++i)
		if( root_config[i] == '1' )
		{
			++N_e;
			K += i;
			if(i<q)
			{
				Kq += i;
				temp1[0] = i;
				temp1++;
			}
		}

	for(auto i=0;i<dk;++i)
	{
		size_t o = 1;
		printf("dk=%d:\n", i);
		edge_counting(N_o-q, N_e-p, K+i-Kq, temp, N_o-q, N_e-p, q, p, o);
		printf("\n");
	}
	return 0;
*/

bool check_root(const Data_basis *temp, const Data_basis N_o, const Data_basis N_p, const Data_basis q, const Data_basis p, size_t& o)
{
    std::vector<char> bin(N_o, 0);
    for(auto i=0;i<N_p;++i)
        bin[temp[i]] = 1;
    int num = 0;
    for(auto i=0;i<q;++i)
    {
        //printf("i = %d, bin = %d\n", i, bin[i]);
        num += (int)bin[i];
    }
    //printf("num = %d\n", num);
    if(num>p) return false;

    for(auto i=q;i<N_o;++i)
    {
        num -= (int)bin[i-q];
        num += (int)bin[i];
        //printf("num = %d\n", num);
        if(num>p) return false;
    }
    printf("%lu:\t", o); 
    for(auto i=0;i<N_o;++i)
        printf("%d ", bin[i]);
    printf("\n");
    ++o;
    return true;
}

int edge_counting(Data_basis N_tot_o, Data_basis N_e, int K, Data_basis *temp,\
                    const Data_basis N_o, const Data_basis N_p, const Data_basis q, const Data_basis p, size_t& o)
{
    if(N_e==0)
    {
        if(K!=0)
            return 0;
        //for(size_t j=0;j<N_p;++j)
        //{
        //    printf("%d ",temp[j]);
        //}
        //printf("\t%lu\n",dict_num(temp,N_p));
        check_root(temp, N_o+q, N_p+p, q, p, o);
        return 0;
    }

    if(K<0)
        return 0;
    for(Data_basis i=N_e-1;i<N_tot_o;++i)
    {
        temp[(size_t)p + N_e - 1] = i + q;
        edge_counting(i, N_e-1, K-i-q, temp, N_o, N_p, q, p, o);
    }
    return 0;
}


