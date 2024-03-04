#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MKL_INT unsigned int
#define MKL_Complex16 complex double
#include<mkl.h>
//#include<mkl_spblas.h>

//--------------------------------------------------------------------------------------------------------------------
//------------------------- declaration ------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

#define ULIMAX 0xFFFFFFFFFFFFFFFF

#include "cQHF_Sphere.h"
#include "Basis/Basis.h"
#include "CGC/CGC.h"
#include "A_list/A_list.h"
#include "Hamiltonian/Hamiltonian.h"

#include "Basis/Basis.c"
#include "CGC/CGC.c"
#include "A_list/A_list.c"
#include "Hamiltonian/Hamiltonian.c"

int main()
{
	FILE *fp = NULL;
	size_t N_e = 4;
	size_t N_phi = 2*N_e;
	int L_z = N_e*(N_e-1)/2;
	printf("init FQHE status: %d\n", init_FQHE(N_phi, N_e));
	printf("create k1 spcae status: %d\n", create_Lz_space(L_z));
	printf("init FQHE status: %d\n", init_FQHE(N_phi, N_e));

	double *pseudo = (double *)calloc(N_e, sizeof(double));
	for(size_t i=0;i<N_e;++i)
		pseudo[i] = (2*i+1)/(4*M_PI);
	double *A_list2 = (double *)calloc((N_e*N_e*N_e*N_e), sizeof(double));
	struct CG2 cg;
	cg.cg = NULL;
	printf("create:\t\t%d\n", createCG2Bdoy(&cg, N_e));
	cal2BodyAlist(&cg, pseudo, A_list2, 1);
	printCG2(&cg);
	return 0;
}

int init_FQHE(size_t N_phi, size_t N_e)
{
	A.N_phi = N_phi;
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
	char *temp =  (char *)malloc((size_t)A.N_e * sizeof(char));
	size_t o = 0;
	st = fermionLzBasisCount(A.N_phi, A.N_e, A.L_z, temp, A.N_e, &o);
	if(st!=0)
		printf("fermionLzBasisCount\t%d\n", st);
	A.basis = (char *)malloc((size_t)o * A.N_e * sizeof(char));
	if(A.basis==NULL)
		return 1;
	A.dim = o;
	o =0;
	st = fermionLzBasis(A.N_phi, A.N_e, A.L_z, temp, A.basis, A.N_e, &o);
	if(st!=0)
		printf("fermionLzBasis\t%d\n", st);
	free(temp);

	A.base_num = (size_t *)malloc(A.dim * sizeof(size_t));
	if(A.base_num==NULL)
		return 1;

	if(A.partition!=NULL)
	{
		free(A.partition);
		A.partition=NULL;
	}

	A.index = (MKL_INT *)malloc(binomial(A.N_phi, A.N_e) * sizeof(MKL_INT));
	memset(A.index, -1, binomial(A.N_phi, A.N_e));

	if(A.index==NULL)
		return 1;

	A.base_num = basis_to_num(A.basis, A.dim, A.N_e, A.index);
	if(A.base_num==NULL)
		return 1;
	char *basis_t = A.basis;
	for(size_t i=0;i<A.dim;++i)
	{
		for(size_t j=0;j<A.N_e;++j)
			printf("%d ",basis_t[j]);
		basis_t += A.N_e;
		printf("\t%u\t%lu\t%u\n",A.base_num[i],i,A.index[A.base_num[i]]);
	}
	return 0;
}