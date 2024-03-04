//--------------------------------------------------------------------------------------------------------------------
//-------------------------<L_z_basis.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

size_t Binomial[50][50] = {0};

size_t binomial(size_t m, size_t n);//C_m^n

size_t dict_num(char *basis, size_t N_e);

size_t gcd(size_t x, size_t y);

int fermionLzBasis(size_t N_phi, size_t N_e, int L_z, char *temp,\
					char *basis, const size_t N_p, size_t *o);

int fermionLzBasisCount(size_t N_phi, size_t N_e, int L_z, char *temp,\
					const size_t N_p, size_t *o);

size_t *basis_to_num(char *basis, size_t dim, size_t N_e, MKL_INT *index);