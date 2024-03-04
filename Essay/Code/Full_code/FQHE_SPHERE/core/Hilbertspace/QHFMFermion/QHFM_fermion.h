
using hash = std::unordered_map<size_t, size_t>;
struct CG2
{
	int u;
	double s;
	double **cg;
	double *head;
};

struct FQHE
{
	size_t N_o;
	size_t N_tot_o;
	size_t N_e;
	size_t L_z;
	int Z_2;
	int PH;

	Data_basis *basis;
	size_t *Z2_num;
	char *Z2_phase;
	char *PH_phase;
	char *Z2PH_length;
	size_t dim;
	size_t Z2dim;
	MKL_INT *PHindex;
	std::unordered_map<size_t, size_t> hashT;
	std::unordered_map<size_t, size_t> PHhashT;

	//struct bipartition* partition; // entanglement partition

	size_t n_data;
	MKL_INT HermitianQ;
	MKL_INT ProjectionQ;

	//struct CG2 cg;
};
struct FQHE A = {
	.N_o = (size_t)0,
	.N_tot_o = (size_t)0,
	.N_e = (size_t)0,
	.L_z = (size_t)0,
	.Z_2 = (int)0,
	.PH = (int)0,
	.basis = NULL,
	.Z2_num = NULL,
	.Z2_phase = NULL,
	.PH_phase = NULL,
	.Z2PH_length = NULL,
	.dim = (size_t)ULIMAX,
	.Z2dim = (size_t)ULIMAX,
	.PHindex = NULL,
	//.partition = NULL,
	.hashT = {},
	.PHhashT = {},
	.n_data = (size_t)ULIMAX,
	.HermitianQ = 0,
	.ProjectionQ = 0,
};
struct FQHE B = A;
/*
	  Z2basis       Z2_phase
	j1 ... jn     {1,-1,2,-2}
		.				.
		.				.
		.				.
		.				.

     PHindex       PH_phase
	   ind        {1,-1,2,-2}
		.				.
		.				.
		.				.
		.				.
*/

struct bipartition
{
	size_t occupy;
	size_t Lz;
	size_t left;
	size_t right;
};

extern "C" int init_FQHE(size_t N_o, size_t N_e);

extern "C" int create_Lz_space(size_t L_z);

extern "C" int create_Lz_Z2_space(size_t L_z, int Z_2);

extern "C" int create_Lz_Z2PH_space(size_t L_z, int Z_2, int PH);// Z_2, PH: -1 or 1

extern "C" void set_Projection(int ProjectionQ);

extern "C" int init_FQHE_type(size_t N_o, size_t N_e, size_t type);
extern "C" int create_Lz_space_type(size_t L_z, size_t type);

extern "C" int get_B_n_o();
extern "C" int get_B_n_e();
extern "C" size_t get_B_Lz_dim();

extern "C" int init_FQHE_B(size_t N_o, size_t N_e);
extern "C" int create_Lz_space_B(size_t L_z);

extern "C" int createA2BdoywithNo(double *pseudo, double *A_list2, int N_o, int num_thread);

extern "C" int createO00AwithNo(double *olist, double *A_list2, int N_o, int num_thread);

extern "C" int createOL0AwithNo(double *olist, double *A_list2, int L, int N_o, int num_thread);

extern "C" int creatennL0AwithNo(double *A_list2, int L, int N_o, int num_thread);

extern "C" int create_nlm_nlm_AwithNo(double *A_list2, double *olist, int L, int N_o, int num_thread);

extern "C" int create_nl1m_nl2m_AwithNo(double *A_list2, int l_1, int l_2, int m, int N_o, int num_thread);

extern "C" int create_delta_AwithNo(double *A_list2, int l, int m, int N_o, int num_thread);

extern "C" int create_deltadelta_ABwithNo(double *A_list2, int l, int m, int la, int lb, int N_o, int num_thread);

extern "C" int create_nlm_AwithNo(double *A_list2, int l, int m, int N_o, int num_thread);

extern "C" int exchangeAB()
{
	struct FQHE C = A;
	A = B;
	B = C;
	return 0;
}