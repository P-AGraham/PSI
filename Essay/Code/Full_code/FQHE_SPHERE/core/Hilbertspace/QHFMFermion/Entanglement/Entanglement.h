//Entanglement.h


//--------------------------------------------------------------------------------------------------------------------
//-------------------------EntanglementMatrix.h-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

int orbitalEntanglementMatrixDimension(MKL_INT *rdmdim, size_t k1, size_t o_left, int num_thread);

int orbitalEntanglementMatrix(MKL_Complex16 *ket_k1, MKL_Complex16 *oem, \
								size_t occupy, size_t ky, size_t dim2, int num_thread);

int orbitalEntanglementSpetrum(MKL_INT *rdmdim, double *s, MKL_Complex16 *ket_k1, \
								size_t s_length, size_t na, int num_thread); //na=o_left
//--------------------------------------------------------------------------------------------------------------------
//-------------------------expectation.h-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
