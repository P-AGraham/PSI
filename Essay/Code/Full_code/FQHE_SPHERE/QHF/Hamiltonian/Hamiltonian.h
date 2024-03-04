

size_t twoBodyCount(MKL_INT *indptr, int num_thread);

size_t twoBodyCount_thread(size_t row);

size_t hopping2Body(unsigned char *basis, unsigned char *bra_i, unsigned char *bra_j, \
					size_t N_phi, size_t N_e, int j1, int j2, int j3, int j4, int *phase);

size_t twoBodyGenerate(double *A_list, MKL_INT *indptr, MKL_INT *indices, \
						double* data, int num_thread);

size_t twoBodyGenerate_thread(size_t row, double *A_list, MKL_INT *indptr, \
								MKL_INT *indices, double *data);