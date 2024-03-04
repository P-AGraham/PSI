#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MKL_Complex16 complex double

#define DataType double
//#define MKL_INT size_t
//#define lapack_int long int
#include <mkl.h>
#include <mkl_spblas.h>
#include <sys/time.h>

//--------------------------------------------------------------------------------------------------------------------
#include "Linear_Algebra/Linalg.h"
#include "SPBLAS/spblas.h"
#include "Lanczos/Lanczos.h"
#include "Lanczos.h"

#include "Linear_Algebra/Linalg.c"
#include "SPBLAS/spblas.c"
#include "Lanczos/Lanczos.c"

//--------------------------------------------------------------------------------------------------------------------

int main()
{
	printf("1 %lu\n", sizeof(MKL_INT));
	FILE *fp = NULL;
	MKL_INT len = 1000000;

	struct timeval start;
	struct timeval end;
	unsigned long timer;

	MKL_INT len_data = 38320; //二进制文件data.bin,indices.bin长度
	MKL_INT len_indptr = 1064+1; //二进制文件indptr.bin长度
	MKL_INT dim = len_indptr - 1;//二进制文件vector.bin,数组out长度

	// 读取indptr.bin并储存在指针indptr指向的内存中
	MKL_INT *indptr = (MKL_INT *)mkl_malloc((size_t)len_indptr * sizeof(MKL_INT), 64);
	//fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/indptr.bin", "rb+");
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/indptr_long.bin", "rb+");
	fread(indptr, sizeof(MKL_INT), len_indptr, fp);
	fclose(fp);
	printf("1\n");

	// 读取indices.bin并储存在指针indices指向的内存中
	MKL_INT *indices = (MKL_INT *) mkl_malloc((size_t)len_data * sizeof(MKL_INT), 64);
	//fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/indices.bin", "rb+");
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/indices_long.bin", "rb+");
	fread(indices, sizeof(MKL_INT), len_data, fp);
	fclose(fp);
	printf("2\n");

	// 读取data.bin并储存在指针data指向的内存中
	DataType *data = (DataType *)mkl_malloc((size_t)len_data * sizeof(DataType), 64);
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/data.bin", "rb+");
	fread(data, sizeof(DataType), len_data, fp);
	fclose(fp);
	printf("3\n");

	// 读取vector.bin并储存在指针vector指向的内存中
	DataType *vector = (DataType *)mkl_malloc((size_t)dim * sizeof(DataType), 64);
	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/vec.bin", "rb+");
	fread(vector, sizeof(DataType), dim, fp);
	fclose(fp);
	printf("4\n");

	// 申请长度为dim的空间储存运算结果并初始化为0
	DataType *out = (DataType *)mkl_calloc((size_t)(dim+1), sizeof(DataType), 64);
	printf("1\n");
	double tol = pow(10, -15);
	create_csr_mkl_d(indptr, indices, data, dim, 0);
	gettimeofday(&start, NULL);
	for(size_t i=0;i<1;++i)
	{
		d_spmv_mkl(vector, out, 1.0, 0.0);
	}
	gettimeofday(&end, NULL);
	printf("time:%lfms\n",(double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_usec-start.tv_usec)/1000);
	destroy_csr_mkl();

	fp = fopen("/home/hu/study/Code/FQHE_Package/FQHE_SPHERE/core/Lanczos/resu.bin", "wb");
	fwrite(out, sizeof(DataType), dim, fp);
	fclose(fp);

	/*gettimeofday(&start, NULL);
	//printf("%d\n",Lanczos_ground(indptr, indices, data, vector, out, dim, tol, 1));
	//Lanczos_quench(indptr, indices, data, vector, out, dim, tol, 1.0, 1);
	gettimeofday(&end, NULL);
	printf("time:%lfms\n",(double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_usec-start.tv_usec)/1000);

	fp = fopen("/home/hu/512/onedriveprivate/study/Code/FQHE_Package/FQHE/core/Lanczos/out_quench.bin", "wb");
	fwrite(out, sizeof(MKL_Complex16), dim, fp);
	fclose(fp);*/

	mkl_free(indptr);
	mkl_free(indices);
	mkl_free(data);
	mkl_free(vector);
	mkl_free(out);

	return 0;
}	

void Lanczos_quench(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
					MKL_Complex16* vector, MKL_Complex16* out, MKL_INT dim, \
					double tol, double dt, MKL_INT HermitianQ)
{
	struct timeval start;
	struct timeval end;
	unsigned long timer;
	gettimeofday(&start, NULL);
	MKL_Complex16	alpha = 1.0;
	MKL_Complex16	beta = 0.0;
	CBLAS_LAYOUT	layout = CblasRowMajor;;
	CBLAS_TRANSPOSE	trans = CblasNoTrans;
	double quench_tol = pow(10,9);

	//Iteration
	if(tol<pow(10,-15))
		tol = pow(10, -15);
	struct Lanczos lan = create_lanczos(indptr, indices, data, dim, tol, vector, HermitianQ);
	for(size_t kk=0;kk<dim;++kk)
	{
		LanczosIteration(&lan);
		quench_tol = quench_tol_Q(&lan, dt);
		printf("tol = %e(%e)\n",quench_tol,lan.tol);
		if(lan.k<0)
		{
			printf("beta = %e\t k = %d\n",lan.beta[lan.n],lan.k);
			break;
		}
		if(quench_tol<tol)
			break;
	}

	// R = [eigv1, eigv2 ...]
	//Cal R exp(-I*dt*\Lambda) R^T (1,0,...,0)^T
	MKL_Complex16* eh = (MKL_Complex16*)mkl_malloc((size_t)lan.n * sizeof(MKL_Complex16), 64);
	for(size_t i=0;i<lan.n;++i)//eh = exp(-I*dt*\Lambda) R^T (1,0,...,0)^T
		eh[i] = cexp(-1.0*I*lan.eigs[i+1]*dt) * lan.eigv[i];

	// V * eh
	MKL_Complex16* psi_t = (MKL_Complex16*)mkl_calloc((size_t) lan.n, sizeof(MKL_Complex16), 64);
	MKL_Complex16* eigv = (MKL_Complex16*)mkl_calloc((size_t)lan.n*lan.n, sizeof(MKL_Complex16), 64);
	for(size_t i=0;i<lan.n*lan.n;++i)
		eigv[i] = (MKL_Complex16)lan.eigv[i];
	cblas_zgemv(layout, trans, lan.n, lan.n, &alpha, eigv, lan.n, eh, 1, &beta, psi_t, 1);
	mkl_free(eh);
	mkl_free(eigv);

	// V = [q1, q2, q3....]
	// V * eh
	trans = CblasTrans;
	cblas_zgemv(layout, trans, lan.n, lan.dim, &alpha, lan.q[1].vec, lan.dim, psi_t, 1, &beta, out, 1);
	
	mkl_free(psi_t);
	destroy_lanczos(&lan);
	gettimeofday(&end, NULL);
	printf("time:%lfms\n",(double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_usec-start.tv_usec)/1000);
	return;
}

int Lanczos_ground(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
					MKL_Complex16* vector, MKL_Complex16* out, MKL_INT dim, \
					double tol, MKL_INT HermitianQ)
{
	struct timeval start;
	struct timeval end;
	unsigned long timer;
	gettimeofday(&start, NULL);
	MKL_Complex16	alpha = 1.0;
	MKL_Complex16	beta = 0.0;
	CBLAS_LAYOUT	layout = CblasRowMajor;;
	CBLAS_TRANSPOSE	trans = CblasNoTrans;
	size_t min_ind, st = -1;

	//Iteration
	if(tol<pow(10,-15))
		tol = pow(10, -15);
	struct Lanczos lan = create_lanczos(indptr, indices, data, dim, tol, vector, HermitianQ);
	for(size_t kk=0;kk<dim;++kk)
	{
		if(kk>500)
			break;
		LanczosIteration(&lan);
		double gtol = first_tol_Q(&lan, &min_ind);
		if(lan.k<0)
		{
			printf("beta = %e\t k = %d\n",lan.beta[lan.n],lan.k);
			// k=-1: beta<tol; k=-2: memmory error;
			st = 1;
			break;
		}
		if(gtol<tol)
		{
			st = 0;
			printf("Lanczos has converged!\n");
			break;
		}
	}
	out[dim] = lan.eigs[min_ind];

	MKL_Complex16 *psi_t = (MKL_Complex16 *)mkl_malloc(lan.n * sizeof(MKL_Complex16), 64);
	for(size_t i=0;i<lan.n;++i)
		psi_t[i] = (MKL_Complex16)lan.eigv[i*lan.n+min_ind-1];

	// V = [q1, q2, q3....]
	// R = [eigv1, eigv2 ...]
	// Cal V * eigv
	trans = CblasTrans;
	cblas_zgemv(layout, trans, lan.n, lan.dim, &alpha, lan.q[1].vec, lan.dim, psi_t, 1, &beta, out, 1);
	mkl_free(psi_t);
	destroy_lanczos(&lan);
	gettimeofday(&end, NULL);
	printf("time:%lfms\n",(double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_usec-start.tv_usec)/1000);
	return st;
}

int Lanczos_ground_tri(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data, \
					MKL_Complex16* vector, MKL_Complex16* out, double *alpha_o, \
					double *beta_o, MKL_INT dim, \
					double tol, MKL_INT HermitianQ)
{
	struct timeval start;
	struct timeval end;
	unsigned long timer;
	gettimeofday(&start, NULL);
	MKL_Complex16	alpha = 1.0;
	MKL_Complex16	beta = 0.0;
	CBLAS_LAYOUT	layout = CblasRowMajor;;
	CBLAS_TRANSPOSE	trans = CblasNoTrans;
	size_t min_ind, st = -1;

	//Iteration
	if(tol<pow(10,-15))
		tol = pow(10, -15);
	struct Lanczos lan = create_lanczos(indptr, indices, data, dim, tol, vector, HermitianQ);
	for(size_t kk=0;kk<dim;++kk)
	{
		if(kk>500)
			break;
		LanczosIteration(&lan);
		double gtol = first_tol_Q(&lan, &min_ind);
		if(lan.k<0)
		{
			printf("beta = %e\t k = %d\n",lan.beta[lan.n],lan.k);
			// k=-1: beta<tol; k=-2: memmory error;
			st = 1;
			break;
		}
		if(gtol<tol)
		{
			st = 0;
			printf("Lanczos has converged!\n");
			break;
		}
	}
	out[dim] = lan.eigs[min_ind];

	MKL_Complex16 *psi_t = (MKL_Complex16 *)mkl_malloc(lan.n * sizeof(MKL_Complex16), 64);
	for(size_t i=0;i<lan.n;++i)
		psi_t[i] = (MKL_Complex16)lan.eigv[i*lan.n+min_ind-1];

	// V = [q1, q2, q3....]
	// R = [eigv1, eigv2 ...]
	// Cal V * eigv
	trans = CblasTrans;
	cblas_zgemv(layout, trans, lan.n, lan.dim, &alpha, lan.q[1].vec, lan.dim, psi_t, 1, &beta, out, 1);
	mkl_free(psi_t);
	
	cblas_dcopy(lan.n, lan.alpha+1, 1, alpha_o, 1);
	cblas_dcopy(lan.n, lan.beta+1, 1, beta_o, 1);
	destroy_lanczos(&lan);
	gettimeofday(&end, NULL);
	printf("time:%lfms\n",(double)(end.tv_sec-start.tv_sec)*1000+(double)(end.tv_usec-start.tv_usec)/1000);
	return st;
}

void create_csr_mkl(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
								 MKL_INT dim, MKL_INT Hermitian)
{
	H = create_sp_csr(indptr, indices, data, dim, Hermitian);return;// H: spblas.h
}

//void spmv_mkl(MKL_Complex16* vec, MKL_Complex16* out, MKL_Complex16 a, MKL_Complex16 b)
void spmv_mkl(MKL_Complex16* vec, MKL_Complex16* out, double ar, double ai, double br, double bi)
{
	//printf("a=%lf+I%lf, b=%lf+I%lf\n", creal(ar+ai*I), cimag(ar+ai*I), creal(br+bi*I), cimag(br+bi*I));
	mkl_sparse_z_mv(H.opp, ar+ai*I, H.mat, H.descr, vec, br+bi*I, out);return;
}

void destroy_csr_mkl()
{
	destroy_sp_csr(H);return;
}








void create_csr_mkl_L2(MKL_INT* indptr, MKL_INT* indices, MKL_Complex16* data,\
								 MKL_INT dim, MKL_INT Hermitian)
{
	L2 = create_sp_csr(indptr, indices, data, dim, Hermitian);return;// H: spblas.h
}

void spmv_mkl_L2(MKL_Complex16* vec, MKL_Complex16* out, double ar, double ai, double br, double bi)
{
	//printf("a=%lf+I%lf, b=%lf+I%lf\n", creal(ar+ai*I), cimag(ar+ai*I), creal(br+bi*I), cimag(br+bi*I));
	mkl_sparse_z_mv(L2.opp, ar+ai*I, L2.mat, L2.descr, vec, br+bi*I, out);return;
}

void destroy_csr_mkl_L2()
{
	destroy_sp_csr(L2);return;
}






void create_csr_mkl_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
								 MKL_INT dim, MKL_INT Hermitian)
{
	printf("create_csr_mkl_d: dim=%ld, Hermitian=%lu\n", dim, Hermitian);

	for(size_t i=0;i<indptr[dim];++i)
	{
		//for(auto j=indptr[i];j<indptr[i+1];++j)
		printf("%lu: %lu, %lf\n", i,indices[i],data[i]);
	}
	H = create_sp_csr_d(indptr, indices, data, dim, Hermitian);return;// H: spblas.h
}

void d_spmv_mkl(double* vec, double* out, double a, double b)
{
	mkl_sparse_d_mv(H.opp, a, H.mat, H.descr, vec, b, out);
	for(size_t i=0;i<3;++i)
		printf("%lf\t%lf\n", vec[i], out[i]);
	return;
}






void create_csr_mkl_L2_d(MKL_INT* indptr, MKL_INT* indices, double* data,\
								 MKL_INT dim, MKL_INT Hermitian)
{
	//printf("create_csr_mkl_L2_d: dim=%lu, Hermitian=%lu\n", dim, Hermitian);
	L2 = create_sp_csr_d(indptr, indices, data, dim, Hermitian);return;// H: spblas.h
}

void d_spmv_mkl_L2(double* vec, double* out, double a, double b)
{
	mkl_sparse_d_mv(L2.opp, a, L2.mat, L2.descr, vec, b, out);return;
}


void my_spmv_z(MKL_Complex16* data, size_t* indices, size_t* indptr, MKL_Complex16* vector, MKL_Complex16* out, double a, double b, size_t dim, size_t n_threads)
{
	size_t i,j;

	# pragma omp parallel for num_threads(n_threads)
	for(i=0;i<dim;i++)
	{	
		out[i] *= b;
		for(j=indptr[i];j<indptr[i+1];j++)
		{
			out[i] += a * data[j] * vector[indices[j]];
		}
	}
}

void my_spmv_d(double* data, size_t* indices, size_t* indptr, double* vector, double* out, double a, double b, size_t dim, size_t n_threads)
{
	size_t i,j;

	# pragma omp parallel for num_threads(n_threads)
	for(i=0;i<dim;i++)
	{	
		out[i] *= b;
		for(j=indptr[i];j<indptr[i+1];j++)
		{
			out[i] += a * data[j] * vector[indices[j]];
		}
	}
}