

struct CG2
{
	int u;
	double s;
	double **cg;
	double *head;
};

struct bipartition
{
	size_t occupy;
	size_t ky;
	size_t left;
	size_t right;
};

struct FQHE
{
	size_t N_phi;
	size_t N_e;
	size_t L_z;

	char *basis;
	size_t *base_num;
	size_t dim;
	MKL_INT *index;

	struct bipartition* partition; // entanglement partition

	size_t n_data;
	MKL_INT HermitianQ;
};
struct FQHE A = {
	.N_phi = (size_t)0,
	.N_e = (size_t)0,
	.L_z = (size_t)0,
	.basis = NULL,
	.base_num = NULL,
	.dim = (size_t)ULIMAX,
	.index = NULL,
	.partition = NULL,
	.n_data = (size_t)ULIMAX,
	.HermitianQ = 0,
};


int init_FQHE(size_t N_phi, size_t N_e);

int create_Lz_space(size_t L_z);