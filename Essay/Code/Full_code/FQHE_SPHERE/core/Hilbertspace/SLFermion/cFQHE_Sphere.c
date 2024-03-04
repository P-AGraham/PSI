#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MKL_INT unsigned int
#define MKL_Complex16 complex double
#define Data_dtype complex double
//#include<mkl.h>
//#include<mkl_spblas.h>

//--------------------------------------------------------------------------------------------------------------------
//------------------------- declaration ------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

#define ULIMAX 0xFFFFFFFFFFFFFFFF

#include "cFQHE_Sphere.h"
#include "Basis/Basis.h"
#include "A_list/A_list.h"
#include "CGC/CGC.h"
#include "Hamiltonian/Hamiltonian.h"

#include "Basis/Basis.c"
#include "A_list/A_list.c"
#include "CGC/CGC.c"
#include "Hamiltonian/Hamiltonian.c"

int main()
{
	FILE *fp = NULL;
	size_t N_o = 19;
	size_t N_e = 7;
	int L_z = 63;
	printf("init FQHE status: %d\n", init_FQHE(N_o, N_e));
	printf("create k1 spcae status: %d\n", create_Lz_space(L_z));

	size_t dk = 17;

	double *pseudo = (double *)calloc(N_o, sizeof(double));
	double *A_list2 = (double *)calloc((N_o*N_o*N_o*N_o), sizeof(double));
	//pseudo[0] = 0.0;
	//pseudo[dk] = sqrt(2.0/(N_phi-1));
	//pseudo[2] = 0.0;
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/C_Test/coulomb19.bin", "rb+");
	fread(pseudo, sizeof(double), 19, fp);
	fclose(fp);


	struct CG2 cg;
	cg.cg = NULL;
	printf("create:\t\t%d\n", createCG2Bdoy(&cg, N_o));
	cal2BodyAlist(&cg, pseudo, A_list2, 1);
	//printCG2(&cg);
	MKL_INT *indptr = (MKL_INT *)calloc(A.dim+1, sizeof(MKL_INT));
	printf("k1_dim = %lu\n",A.dim);
	printf("pre-work for parallel; nnz = : %lu\n", twoBodyCount(indptr, 6));

	MKL_INT *indices = (MKL_INT *)malloc(indptr[A.dim] * sizeof(MKL_INT));
	Data_dtype *data = (Data_dtype *)malloc(indptr[A.dim] * sizeof(Data_dtype));
	printf("hamiltonian_generate status: %lu\n", twoBodyGenerate(A_list2, indptr, indices, data, 6));

	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/C_Test/Ne7/indptr17.bin", "wb");
	fwrite(indptr, sizeof(MKL_INT), A.dim+1, fp);
	fclose(fp);
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/C_Test/Ne7/indices17.bin", "wb");
	fwrite(indices, sizeof(MKL_INT), indptr[A.dim], fp);
	fclose(fp);
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/C_Test/Ne7/data17.bin", "wb");
	fwrite(data, sizeof(double), indptr[A.dim], fp);
	fclose(fp);
	printf("init FQHE status: %d\n", init_FQHE(N_o, N_e));
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/C_Test/Ne7/A_list15.bin", "wb");
	fwrite(A_list2, sizeof(double), (N_o*N_o*N_o*N_o), fp);
	fclose(fp);

	free(data);
	free(indices);
	free(indptr);
	printf("destroy:\t%d\n", destoryCG2(&cg));
	free(pseudo);
	free(A_list2);
	return 0;
}

int init_FQHE(size_t N_o, size_t N_e)
{
	A.N_o = N_o;
	A.N_e = N_e;

	A.L_z = 0;
	A.n_data = ULIMAX;

	if(A.basis!=NULL)
	{
		free(A.basis);
		A.basis = NULL;
	}

	if(A.base_num!=NULL)
	{
		free(A.base_num);
		A.base_num = NULL;
	}

	A.dim = ULIMAX;

	if(A.index!=NULL)
	{
		free(A.index);
		A.index = NULL;
	}

	if(A.partition!=NULL)
	{	
		free(A.partition);
		A.partition=NULL;
	}
	return 0;
}

int createA2Bdoy(double *pseudo, double *A_list2, int num_thread)
{
	struct CG2 cg;
	cg.cg = NULL;
	int st = 0;
	st = createCG2Bdoy(&cg, A.N_o);
	if(st!=0)
	{
		printf("createCG2Bdoy\t%d\n", st);
		return st;
	}
	cal2BodyAlist(&cg, pseudo, A_list2, num_thread);
	destoryCG2(&cg);
	return st;
}


int createA2BdoywithNo(double *pseudo, double *A_list2, int N_o, int num_thread)
{
	struct CG2 cg;
	cg.cg = NULL;
	int st = 0;
	st = createCG2Bdoy(&cg, N_o);
	if(st!=0)
	{
		printf("createCG2Bdoy\t%d\n", st);
		return st;
	}
	cal2BodyAlist(&cg, pseudo, A_list2, num_thread);
	destoryCG2(&cg);
	return st;
}

int create_Lz_space(size_t L_z)
{
	if(A.L_z == L_z)
		return 0;
	A.L_z = L_z;

	if(A.basis!=NULL)
	{
		free(A.basis);
		A.basis = NULL;
	}
	if(A.base_num!=NULL)
	{
		free(A.base_num);
		A.base_num = NULL;
	}
	if(A.index!=NULL)
	{
		free(A.index);
		A.index = NULL;
	}

	int st;
	size_t o = 0;
	char *temp =  (char *)malloc((size_t)A.N_e * sizeof(char));
	st = fermionLzBasisCount(A.N_o, A.N_e, A.L_z, temp, A.N_e, &o);
	if(st!=0)
	{
		printf("fermionLzBasisCount\t%d\n", st);
		return st;
	}
	A.basis = (char *)malloc((size_t)o * A.N_e * sizeof(char));
	if(A.basis==NULL)
		return 1;
	A.dim = o;
	o =0;
	st = fermionLzBasis(A.N_o, A.N_e, A.L_z, temp, A.basis, A.N_e, &o);
	if(st!=0)
	{
		printf("fermionLzBasis\t%d\n", st);
		return st;
	}
	free(temp);

	A.base_num = (size_t *)malloc(A.dim * sizeof(size_t));
	if(A.base_num==NULL)
		return 1;

	if(A.partition!=NULL)
	{
		free(A.partition);
		A.partition=NULL;
	}

	A.index = (MKL_INT *)malloc(binomial(A.N_o, A.N_e) * sizeof(MKL_INT));
	memset(A.index, -1, binomial(A.N_o, A.N_e));

	if(A.index==NULL)
		return 1;

	A.base_num = basis_to_num(A.basis, A.dim, A.N_e, A.index);
	if(A.base_num==NULL)
		return 1;
	/*char *basis_t = A.basis;
	for(size_t i=0;i<A.dim;++i)
	{
		for(size_t j=0;j<A.N_e;++j)
			printf("%d ",basis_t[j]);
		basis_t += A.N_e;
		printf("\t%u\t%lu\t%u\n",A.base_num[i],i,A.index[A.base_num[i]]);
	}*/
	return 0;
}


void copy_basis(char* basis)
{
	for(size_t i=0;i<A.dim*A.N_e;++i)
		basis[i] = A.basis[i];
	return;
}

int get_n_o(){return (int)A.N_o;}
int get_n_e(){return (int)A.N_e;}
void set_Hermitian(MKL_INT HermitianQ){A.HermitianQ=HermitianQ;}
MKL_INT get_Hermitian(){return A.HermitianQ;}
size_t get_Lz_dim(){return A.dim;}