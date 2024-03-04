
#include "quicksort.h"

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<twobody.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" size_t twoBodyCount(MKL_INT *indptr, int num_thread);

size_t twoBodyCount_thread(Data_basis *k1_basis, size_t k1_dim, \
					size_t N_o, size_t N_e, size_t i);

double cdw_energy(Data_basis *basis, double *cdw, size_t N_phi, size_t N_e);

size_t tight_bond_count(Data_basis *basis, size_t N_phi, size_t N_e);

size_t tight_up(size_t j1, size_t N_phi);

size_t tight_down(size_t j1, size_t N_phi);

size_t layerQ(size_t j1, size_t j2, size_t N_phi);

extern "C" size_t twoBodyGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
					double *A1_list, double *A12_list, double *cdw, double t, int num_thread);

size_t twoBodyGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
					double *A1_list, double *A12_list, double *cdw, \
					Data_basis *k1_basis, const hash &hashT, double t, \
					size_t k1_dim, size_t N_o, size_t N_e, size_t row);

size_t tight_bond(Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
				const hash & hashT, MKL_INT *indices_thread, Data_dtype *data_thread, \
				double t, size_t N_phi, size_t N_e, size_t count);

size_t tight_bond_hopping(Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
							size_t N_e, size_t N_phi, size_t j1, size_t j2, int *phase);

size_t intra_hopping(const Data_basis *basis, Data_basis *bra_i, Data_basis *bra_j, \
					size_t N_phi, size_t N_e, int j1, int j2, int j3, int j4, int *phase);
size_t inter_hopping(Data_basis *basis, Data_basis *bra_i, Data_basis *bra_j, \
					size_t N_orbit, size_t N_e, int j1, int j2, int j3, int j4, int *phase);

size_t intra_hopping(const Data_basis *basis, Data_basis *bra_i, Data_basis *bra_j, \
					size_t N_orbit, size_t N_e, int j1, int j2, int j3, int j4, int *phase)
{
	// j4 > j3, j1 > j2

	// initialization
	memset(bra_i, 0, (size_t)N_orbit * sizeof(Data_basis));

	// translation occupy -> binary
	for(size_t i=0;i<N_e;++i)
		++bra_i[basis[i]];
	*phase = 2;

	// annihilation
	for(size_t i=j3+1;i<j4;++i)
		*phase += bra_i[i];
	--bra_i[j4];
	--bra_i[j3];

	// creation
	for(size_t i=j2+1;i<j1;++i)
		*phase += bra_i[i];
	++bra_i[j2];
	++bra_i[j1];

	*phase = (*phase&1)?-1:1;
	
	// translation binary -> occupy
	size_t k = 0;
	for(size_t i=0;i<N_orbit;++i)
	{
		if(bra_i[i])
		{
			bra_j[k] = i;
			++k;
		}
	}

	return dict_num(bra_j, N_e);
}


//--------------------------------------------------------------------------------------------------------------------
//-------------------------<L2.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

bool Lplus(Data_basis *bra_i, size_t j, size_t N_o);// L_j^+

bool Lminus(Data_basis *bra_i, size_t j, size_t N_o);// L_j^-

size_t Lplus(const Data_basis *bra_i, Data_basis *bra_j, size_t j, size_t N_o, size_t N_e);// L_j^+

size_t Lminus(const Data_basis *bra_i, Data_basis *bra_j, size_t j, size_t N_o, size_t N_e);// L_j^-

bool LplusQ(const Data_basis *bra_i, size_t j, size_t N_o);// L_j^+

bool LminusQ(const Data_basis *bra_i, size_t j, size_t N_o);// L_j^-

extern "C" size_t L2Count(MKL_INT *indptr, int num_thread);

size_t L2Count_thread(Data_basis *k1_basis, size_t k1_dim, \
                    size_t N_o, size_t N_e, size_t i);

extern "C" size_t L2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);

size_t L2Generate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
						Data_basis *k1_basis, const hash &hashT, \
						size_t k1_dim, size_t N_o, size_t N_e, size_t row);



//--------------------------------------------------------------------------------------------------------------------
//-------------------------<twobody_Z2.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" size_t twoBodyZ2Count(MKL_INT *indptr, int num_thread);

size_t twoBodyZ2Count_thread(Data_basis *k1_basis, const hash & hashT, char *Z2_phase, int Z_2, size_t k1_dim, \
                    size_t N_o, size_t N_e, size_t i);

size_t tight_bond_Z2_count(const Data_basis *basis, Data_basis *bra_i, Data_basis *bra_j, const hash & hashT, size_t N_o, size_t N_e);


extern "C" size_t twoBodyZ2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double t, int num_thread);

size_t twoBodyZ2Generate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, \
                    Data_basis *k1_basis, const hash & hashT, char* Z2_phase, int Z_2, double t, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row);


size_t tight_bond_Z2(const Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
                const hash & hashT, const char* Z2_phase, int Z_2, MKL_INT *indices_thread, Data_dtype *data_thread, \
                double t, size_t N_phi, size_t N_e, size_t count);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<L2_Z2.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" size_t L2Z2Count(MKL_INT *indptr, int num_thread);

size_t L2Z2Count_thread(Data_basis *k1_basis, const hash & hashT, char* Z2_phase, int Z_2, size_t k1_dim, \
                    size_t N_o, size_t N_e, size_t i);

extern "C" size_t L2Z2Generate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);

size_t L2Z2Generate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, Data_basis *k1_basis, \
                            const hash & hashT, char* Z2_phase, int Z_2, size_t k1_dim, size_t N_o, size_t N_e, size_t row);


//--------------------------------------------------------------------------------------------------------------------
//-------------------------<twobody_Z2PH.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------


extern "C" size_t twoBodyZ2PHCount(MKL_INT *indptr, int num_thread);

size_t twoBodyZ2PHCount_thread(const Data_basis *k1_basis, const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                             const char* Z2_phase, const char* Z2PH_phase, const char* PH_length, int Z_2, int PH,\
                              size_t k1_dim, size_t N_o, size_t N_e, size_t i);

size_t tight_bond_Z2PH_count(const Data_basis *basis, Data_basis *bra_i, Data_basis *bra_j, const size_t* PHindex, const hash & hashT,\
                        const hash & PHhashT, const char* Z2_phase, const char* PH_phase,  const char* Z2PH_length, size_t N_o, size_t N_e);

extern "C" size_t twoBodyZ2PHGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double t, int sortQ, int num_thread);

size_t twoBodyZ2PHGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, double *A1_list, double *A12_list, \
                    Data_basis *k1_basis, const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                    const char* Z2_phase, const char* PH_phase, const char* Z2PH_length, int Z_2, int PH, double t, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row);

size_t tight_bond_Z2PH(const Data_basis *basis_t, Data_basis *bra_i, Data_basis *bra_j, \
                 const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                    const char* Z2_phase, const char* PH_phase, const char* Z2PH_length, int Z_2, int PH, \
                    MKL_INT *indices_thread, Data_dtype *data_thread, double t, size_t N_o, size_t N_e, size_t row, int ketlength);

//--------------------------------------------------------------------------------------------------------------------
//-------------------------<L2_Z2PH.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" size_t L2Z2PHCount(MKL_INT *indptr, int num_thread);

size_t L2Z2PHCount_thread(Data_basis *k1_basis, const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                             const char* Z2_phase, const char* PH_phase, const char* Z2PH_length, int Z_2, int PH, \
                             size_t k1_dim, size_t N_o, size_t N_e, size_t i);

extern "C" size_t L2Z2PHGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);

size_t L2Z2PHGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    Data_basis *k1_basis, const size_t* PHindex, const hash & hashT, const hash & PHhashT,\
                    const char* Z2_phase, const char* PH_phase, const char* Z2PH_length, int Z_2, int PH, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row);


//--------------------------------------------------------------------------------------------------------------------
//-------------------------<nlm.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

bool hopping_UptoDown(Data_basis* bin_basis, Data_basis* resu_basis, size_t m, size_t N_e, size_t N_o, char& phase);
bool hopping_DowntoUp(Data_basis* bin_basis, Data_basis* resu_basis, size_t m, size_t N_e, size_t N_o, char& phase);
bool hopping_UptoDown(Data_basis* bin_basis, Data_basis* resu_basis, int m1, int m, size_t N_e, size_t N_o, char& phase);
bool hopping_DowntoUp(Data_basis* bin_basis, Data_basis* resu_basis, int m1, int m, size_t N_e, size_t N_o, char& phase);
/*
    layerm1      layerm    hopping_type            hopping                        hoppping_in_code
        0           0       Up_to_Up       c^dagger_{up, m1}c_{up,m1-m}        c^dagger_{m1}c_{m1-m}
        1           0       Up_to_Dn       c^dagger_{dn, m1}c_{up,m1-m}      c^dagger_{m1+N_o}c_{m1-m}
        0           1       Dn_to_Up       c^dagger_{up, m1}c_{dn,m1-m}      c^dagger_{m1}c_{m1-m+N_o}
        1           1       Dn_to_Dn       c^dagger_{dn, m1}c_{dn,m1-m}     c^dagger_{m1+N_o}c_{m1-m+N_o}
*/
bool hopping_LayertoLayer(Data_basis* bin_basis, Data_basis* resu_basis, int m1, int m, int layerm1, int layerm,
                             size_t N_e, size_t N_o, char& phase);
void restore_LayertoLayer(Data_basis* bin_basis, int m1, int m, int layerm1, int layerm, size_t N_o);

extern "C" size_t n00ACount(MKL_INT *indptr, Data_dtype *A_matirx, int num_thread);

size_t n00ACount_thread(const Data_dtype *A_matirx, Data_basis *basis, const hash & hashT, size_t k1_dim, size_t N_o, size_t N_e, size_t i);

extern "C" size_t n00AGenerate(Data_dtype *A_matirx, MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);

size_t n00AGenerate_thread(const Data_dtype *A_matirx, MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                            Data_basis *basis, const hash & hashT, size_t k1_dim, size_t N_o, size_t N_e, size_t row);

extern "C" size_t D00ACount(MKL_INT *indptr, Data_dtype *A_matirx, Data_dtype *B_matirx, int num_thread);
size_t D00ACount_thread(const Data_dtype *A_matirx, const Data_dtype *B_matirx, Data_basis *basis, const hash & hashT,
                         size_t k1_dim, size_t N_o, size_t N_e, size_t i);

extern "C" size_t D00AGenerate(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, 
											MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, int num_thread);
size_t D00AGenerate_thread(const Data_dtype *A_matirx, const Data_dtype *B_matirx, const Data_dtype *A_list, 
                            MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, Data_basis *basis, const hash & hashT, 
                            size_t k1_dim, size_t N_o, size_t N_e, size_t row);

extern "C" void D00AOp(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, int num_thread);
void D00AOpGenerate_thread(const Data_dtype *A_matirx, const Data_dtype *B_matirx, const Data_dtype *A_list, 
                            const Data_dtype *vec, Data_dtype *resu, Data_basis *basis, const hash & hashT, 
                            size_t k1_dim, size_t N_o, size_t N_e, size_t row);

extern "C" void nl0AOp(Data_dtype *A_matirx, Data_dtype *vec, Data_dtype *resu, int l, int num_thread);
void nl0AOpGenerate_thread(const Data_dtype *A_matirx, int l, const struct CG2 &cg, const Data_dtype *vec, Data_dtype *resu, Data_basis *basis, const hash & hashT, 
                            size_t k1_dim, size_t N_o, size_t N_e, size_t row);


void my_spmv_d(double* data, size_t* indices, size_t* indptr, double* vector, double* out, double a, double b, size_t dim, size_t n_threads);
void my_spmv_d(double* data, size_t* indices, size_t* indptr, double* vector, double* out, size_t dim, size_t n_threads);



//--------------------------------------------------------------------------------------------------------------------
//-------------------------<General_operator.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

long creation_two(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e, int j1, int j2, int& phase);
long destroy_two(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e, int j1, int j2, int& phase);

// daggerQ: 
//          0: a_{m-m_1}a_{m_1}  
//          1: a_{m_1}^daggera_{m-m_1}^dagger
extern "C" void pairingOp(Data_dtype *A_list, Data_dtype *A_matirx, Data_dtype *vec, Data_dtype *resu, int l, int m, 
                            int daggerQ, int num_thread);

//a_{m-m_1}a_{m_1} 
void pairingOpGenerate_thread(const Data_dtype *A_list, const Data_dtype *A_matirx, const Data_dtype *vec, Data_dtype *resu, int m,
                            Data_basis *rbasis, const hash & rhashT, size_t rk1_dim, size_t rN_o, size_t rN_e, 
                            Data_basis *lbasis, const hash & lhashT, size_t lk1_dim, size_t lN_o, size_t lN_e, int daggerQ, size_t row);

Data_dtype checkLayer(const Data_dtype *A_matirx, int j1, int j2, int N_o, int layers, bool daggerQ);

long int pairing_hopping(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e,
                         int j1, int j2, int j3, int j4, int& phase);

extern "C" void pairingABlmOp(Data_dtype *A_matirx, Data_dtype *B_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, 
                                int l, int m, int la, int lb, int num_thread);

void pairingABlmOpGenerate_thread(const Data_dtype *A_list, const Data_dtype *A_matirx,  const Data_dtype *B_matirx, 
                            const Data_dtype *vec, Data_dtype *resu, int l, int m, int la, int lb,
                            Data_basis *rbasis, const hash & rhashT, size_t rk1_dim, size_t N_o, size_t N_e, 
                            Data_basis *lbasis, const hash & lhashT, size_t lk1_dim, size_t row);

extern "C" void twoOpGeneral(Data_dtype *A_matirx, Data_dtype *A_list, Data_dtype *vec, Data_dtype *resu, int m, int num_thread);

void twoOpGeneralGenerate_thread(const Data_dtype *A_list, const Data_dtype *A_matirx, const Data_dtype *vec, Data_dtype *resu, int m,
                            Data_basis *rbasis, const hash & rhashT, size_t rk1_dim, size_t N_o, size_t N_e, 
                            Data_basis *lbasis, const hash & lhashT, size_t lk1_dim, size_t row);

long int nlm_hopping(const Data_basis* bin_basis, Data_basis* bra, Data_basis* ket, int N_o, int N_e,
                         int j1, int j2, int& phase);



//--------------------------------------------------------------------------------------------------------------------
//-------------------------<lineDefect.c>-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

extern "C" size_t defectLineGenerate(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, double *A1_list, double *A12_list, double *cdw, double* defect_list, double t, int num_thread);

size_t defectLineGenerate_thread(MKL_INT *indptr, MKL_INT *indices, Data_dtype *data, \
                    double *A1_list, double *A12_list, double *cdw, double* defect_list,\
                    Data_basis *k1_basis, const hash & hashT, double t, \
                    size_t k1_dim, size_t N_o, size_t N_e, size_t row);