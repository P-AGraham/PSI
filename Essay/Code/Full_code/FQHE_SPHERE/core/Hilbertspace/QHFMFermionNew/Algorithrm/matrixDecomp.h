
#pragma once

#include<iostream>
#include "matrixWrapper.h"


namespace mkl{
#include<mkl.h>

//----------------------------------------------------------
//--------------------------svd-----------------------------
//----------------------------------------------------------

long dgesvd(double* matrix, double* s, long m, long n)
{
	char jobu = 'N';//the contents of matrix are destroyed!!
	char jobv = 'N';//the contents of matrix are destroyed!!
	double *u = NULL;
	double *vt= NULL;
	int layout = LAPACK_ROW_MAJOR;
	double *superb = (double *)mkl_calloc(m*n, sizeof(double), 64);

	long lda = n;
	long ldu = m;
	long ldvt = n;

	lapack_int info = LAPACKE_dgesvd(layout, jobu, jobv, m, n, matrix, lda, s, u, ldu, vt, ldvt, superb);
	if(info!=0)
		printf("Unconverged! Info=%lld\n", info);

	mkl_free(superb);
	return info;
}


// If jobz≠'O', the contents of matrix are destroyed!!
long dgesdd(double* matrix, double* s, long m, long n)
{
	char jobz = 'N';
	double *u = NULL;
	double *vt= NULL;
	int layout = LAPACK_ROW_MAJOR;

	long lda = n;
	long ldu = m;
	long ldvt = n;

	lapack_int info = LAPACKE_dgesdd(layout, jobz, m, n, matrix, lda, s, u, ldu, vt, ldvt);
	if(info!=0)
		printf("dgesdd unconverged! Info=%lld\n", info);
	
	return info;
}

//----------------------------------------------------------
//-------------------------eigen----------------------------
//----------------------------------------------------------

long dsyevd(matrixWrapper& A, double* w, bool vecQ, bool upQ)
{
	if(A.dimRow_!=A.dimCol_)
	{
        printf("dsyevd: Mismatch of matrix dimenisons: (%ld, %ld)\n", A.dimRow_, A.dimCol_);
        throw std::runtime_error("dsyevd: Mismatch of matrix dimenisons\n");
	}
	auto matrix_layout = LAPACK_ROW_MAJOR;
	char jobz = vecQ?'V':'N';
	char uplo =  upQ?'U':'L';
	lapack_int n = A.dimRow_;
	lapack_int lda = n;
	double* a = A.matrix_;

	auto info = LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);       
	if(info!=0)
		printf("dsyevd unconverged! Info=%lld\n", info);
	return info;
}
long dsyevd(matrixWrapper& A, double* w)
{
	return dsyevd(A, w, true, true);
}

//lapack_int LAPACKE_dsyev (int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);
long dsyev(matrixWrapper& A, double* w, bool vecQ, bool upQ)
{
	if(A.dimRow_!=A.dimCol_)
	{
        printf("dsyev: Mismatch of matrix dimenisons: (%ld, %ld)\n", A.dimRow_, A.dimCol_);
        throw std::runtime_error("dsyev: Mismatch of matrix dimenisons\n");
	}
	auto matrix_layout = LAPACK_ROW_MAJOR;
	char jobz = vecQ?'V':'N';
	char uplo =  upQ?'U':'L';
	lapack_int n = A.dimRow_;
	lapack_int lda = n;
	double* a = A.matrix_;

	auto info = LAPACKE_dsyev(matrix_layout, jobz, uplo, n, a, lda, w);       
	if(info!=0)
		printf("dsyevd unconverged! Info=%lld\n", info);
	return info;
}
long dsyev(matrixWrapper& A, double* w)
{
	return dsyev(A, w, true, true);
}

}