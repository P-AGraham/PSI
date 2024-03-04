
#pragma once

#include "matrixWrapper.h"

namespace mkl{
#include<mkl.h>

// C = alpha*op(A)*op(B) + beta*C
void matrixMult(const matrixWrapper& A, const matrixWrapper& B, matrixWrapper& C, bool transA, bool transB, double alpha, double beta)
{
    long dimARow, dimBRow, dimCRow, dimACol, dimBCol, dimCCol;

    dimARow = (!transA)?A.dimRow_:A.dimCol_;
    dimACol = ( transA)?A.dimRow_:A.dimCol_;

    dimBRow = (!transB)?B.dimRow_:B.dimCol_;
    dimBCol = ( transB)?B.dimRow_:B.dimCol_;

    dimCRow = C.dimRow_;
    dimCCol = C.dimCol_;

    if(dimACol!=dimBRow || dimARow!=dimCRow || dimBCol!=dimCCol)
    {
        printf("Mismatch of matrix dimenisons: (%ld, %ld)*(%ld, %ld) != (%ld, %ld)\n", dimARow, dimACol, dimBRow, dimBCol, dimCRow, dimCCol );
        throw std::runtime_error("Mismatch of matrix dimenisons\n");
    }

    const double* a = A.matrix_;
    const double* b = B.matrix_;
    double* c = C.matrix_;

    auto Layout = CblasRowMajor;
    auto transa = (transA)?CblasTrans:CblasNoTrans;
    auto transb = (transB)?CblasTrans:CblasNoTrans;

    // (m, k)*(k, n) = (m, n)
    long m = dimARow;
    long k = dimACol;
    long n = dimBCol;

    MKL_INT lda = std::max( ((transA)?m:k), (long)1);
    MKL_INT ldb = std::max( ((transB)?k:n), (long)1);
    MKL_INT ldc = std::max( n, (long)1);

    //             1       2       3    4  5  6    7    8   9  10  11    12  13  14
    cblas_dgemm(Layout, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

    return;
}

// C = op(A)*op(B);
void matrixMult(const matrixWrapper& A, const matrixWrapper& B, matrixWrapper& C, bool transA, bool transB)
{
    matrixMult(A, B, C, transA, transB, 1.0, 0.0);
}

};