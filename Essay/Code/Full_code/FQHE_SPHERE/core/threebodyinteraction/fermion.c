#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

void cal_lambda(double *lambda, double *qx, double *qy, int bound, size_t ax, double g11, double g12, double g22);
void A_Three_Body_0_summary(complex double *A_list, complex double *A_list_0, int N_phi);
inline size_t A_summary(int N_phi, int j1, int j2, int j3, int j4, int j5, int j6);
void A_Three_Body_0_thread(complex double *A_list, double *lambda, const complex double *ftc,\
							int bound, int N_phi, int j16);

void A_Three_Body_0(complex double *A_list, complex double *A_list_0, double *qx, double *qy, \
					int bound, int N_phi, int num_thread, double g11, double g12, double g22);

inline size_t A_summary_0(int N_phi, int j1, int j2, int j3, int j4, int j5, int j6);


void A_Three_Body_0_thread(complex double *A_list, double *lambda, const complex double *ftc,\
							int bound, int N_phi, int j16)
{
	size_t stride, s, sb, ya, yb, l1, l2, l3, l4;
	int temp;
	s = 2*N_phi - 1;
	stride = s * s;
	sb = 2*bound + 1;
	l3 = sb;
	l2 = sb * sb;
	l1 = sb * sb * sb;
	double rho, q1, q2, q12;
	complex double beta;

	for(int j43=N_phi-1;j43<s;++j43)
		for(int j62=N_phi-1;j62<s;++j62)
			for(int j24=N_phi-1;j24<s;++j24)
			{
				beta = 0;
				for(int ax=(1+j16)%N_phi;ax<sb;ax+=N_phi)
				{
					for(int bx=(1+j43)%N_phi;bx<sb;bx+=N_phi)
					{
						for(int ay=0;ay<sb;++ay)
						{
							for(int by=0;by<sb;++by)
							{
								l4 = (size_t)ax*l1 + (size_t)bx*l2 + (size_t)ay*l3 + by;
								temp = ((ay-bound)*(2*(j62-N_phi+1)+bx-bound) + (by-bound)*(2*(j24-N_phi+1)+ax-bound))%(2*N_phi);
								temp += (temp<0)?(2*N_phi):0;
								beta += /*cexp(temp)*/ftc[temp] * lambda[l4];
							}
						}
					}
				}
				A_list[j43*stride+j62*s+j24] = beta;
			}

	return;
}


void A_Three_Body_0(complex double *A_list, complex double *A_list_0, double *qx, double *qy, \
					int bound, int N_phi, int num_thread, double g11, double g12, double g22)
{
	size_t stride, s, d;
	d = 2*bound + 1;
	s = 2*N_phi - 1;
	stride = s * s * s;
	complex double *ftc = (complex double*)malloc((size_t) 2 * N_phi * sizeof(complex double));

	for(size_t i=0;i<2*N_phi;++i)
		ftc[i] = cexp(-1.0*I*M_PI*i/N_phi);

	double *lambda = (double *)malloc(d*d*d*d * sizeof(double));
	# pragma omp parallel for num_threads(num_thread)
	for(size_t ax=0;ax<d;++ax)
		cal_lambda(lambda+ax*d*d*d, qx, qy, bound, ax, g11, g12, g22);
	
	num_thread = (num_thread>N_phi)?N_phi:num_thread;

	# pragma omp parallel for num_threads(num_thread)
	for(int i=N_phi-1;i<s;++i)
		A_Three_Body_0_thread(A_list+i*stride, lambda, ftc, bound, N_phi, i);

	free(lambda);
	free(ftc);

	A_Three_Body_0_summary(A_list, A_list_0, N_phi);

	return;
}

void cal_lambda(double *lambda, double *qx, double *qy, int bound, size_t ax, double g11, double g12, double g22)
{
	size_t stride, d;
	d = 2*bound + 1;
	stride = d * d;
	double qax, qay, qbx, qby, q1, q2, q12, q1o, q2o, q12o;
	qax = qx[ax];

	for(size_t bx=0;bx<d;++bx)
	{
		qbx = qx[bx];
		for(size_t ay=0;ay<d;++ay)
		{
			qay = qy[ax*d+ay];
			for(size_t by=0;by<d;++by)
			{
				qby = qy[bx*d+by];
				q1 = g11*qax*qax + 2.0*g12*qax*qay + g22*qay*qay;
				q2 = g11*qbx*qbx + 2.0*g12*qbx*qby + g22*qby*qby;
				q12=  g11*(qax-qbx)*(qax-qbx) + 2.0*g12*(qax-qbx)*(qay-qby) + g22*(qay-qby)*(qay-qby);
				q1o = qax*qax + qay*qay;
				q2o = qbx*qbx + qby*qby;
				q12o= (qax-qbx)*(qax-qbx) + (qay-qby)*(qay-qby);
				lambda[bx*stride+ay*d+by] = q1o * q1o * q12o * exp(-(q1+q2+q12)/4.0);
				//printf("%lf\t",lambda[bx*stride+ay*d+by]);
			}
		}
	}
	return;
}

void A_Three_Body_0_summary(complex double *A_list, complex double *A_list_0, int N_phi)
{
	size_t stride, s, d, ad;
	s = 2*N_phi - 1;
	d = s * s;
	stride = s * s * s;
	int j1, j2, j3, j4, j5, j6, J123, J456, J45, J32;

	for(size_t u=0;u<N_phi;++u)
	{
		j4 = u;
		for(size_t v=0;v<N_phi;++v)
		{
			j5 = v;
			J45 = j4 + j5;
			if(J45>=N_phi)
				J45 -= N_phi;
			for(size_t w=0;w<N_phi;++w)
			{
				j6 = w;
				J456 = J45 + j6;
				if(J456>=N_phi)
					J456 -= N_phi;
				for(size_t r=0;r<N_phi;++r)
				{
					j3 = r;
					for(size_t s=0;s<N_phi;++s)
					{
						j2 = s;
						J32 = j2 + j3;
						if(J32>=N_phi)
							J32 -= N_phi;
						for(size_t t=0;t<N_phi;++t)
						{
							j1 = t;
							J123 = J32 + j1;
							if(J123>=N_phi)
								J123 -= N_phi;
							if(J123!=J456)
								continue;
							ad = A_summary_0(N_phi, j1, j2, j3, j4, j5, j6);
							if(A_list_0[ad]!=0.0)
								continue;

							A_list_0[ad] = A_list[A_summary(N_phi, j1, j2, j3, j4, j5, j6)];
							//printf("%d %d %d %d %d %d %d %d\n",j1,j2,j3,j4,j5,j6,J123,J456);
							A_list_0[ad] -= A_list[A_summary(N_phi, j1, j3, j2, j4, j5, j6)];
							A_list_0[ad] += A_list[A_summary(N_phi, j2, j3, j1, j4, j5, j6)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j2, j1, j3, j4, j5, j6)];
							A_list_0[ad] += A_list[A_summary(N_phi, j3, j1, j2, j4, j5, j6)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j3, j2, j1, j4, j5, j6)];

							A_list_0[ad] -= A_list[A_summary(N_phi, j1, j2, j3, j4, j6, j5)];
							A_list_0[ad] += A_list[A_summary(N_phi, j1, j3, j2, j4, j6, j5)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j2, j3, j1, j4, j6, j5)];
							A_list_0[ad] += A_list[A_summary(N_phi, j2, j1, j3, j4, j6, j5)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j3, j1, j2, j4, j6, j5)];
							A_list_0[ad] += A_list[A_summary(N_phi, j3, j2, j1, j4, j6, j5)];

							A_list_0[ad] += A_list[A_summary(N_phi, j1, j2, j3, j5, j6, j4)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j1, j3, j2, j5, j6, j4)];
							A_list_0[ad] += A_list[A_summary(N_phi, j2, j3, j1, j5, j6, j4)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j2, j1, j3, j5, j6, j4)];
							A_list_0[ad] += A_list[A_summary(N_phi, j3, j1, j2, j5, j6, j4)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j3, j2, j1, j5, j6, j4)];

							A_list_0[ad] -= A_list[A_summary(N_phi, j1, j2, j3, j5, j4, j6)];
							A_list_0[ad] += A_list[A_summary(N_phi, j1, j3, j2, j5, j4, j6)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j2, j3, j1, j5, j4, j6)];
							A_list_0[ad] += A_list[A_summary(N_phi, j2, j1, j3, j5, j4, j6)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j3, j1, j2, j5, j4, j6)];
							A_list_0[ad] += A_list[A_summary(N_phi, j3, j2, j1, j5, j4, j6)];

							A_list_0[ad] += A_list[A_summary(N_phi, j1, j2, j3, j6, j4, j5)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j1, j3, j2, j6, j4, j5)];
							A_list_0[ad] += A_list[A_summary(N_phi, j2, j3, j1, j6, j4, j5)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j2, j1, j3, j6, j4, j5)];
							A_list_0[ad] += A_list[A_summary(N_phi, j3, j1, j2, j6, j4, j5)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j3, j2, j1, j6, j4, j5)];

							A_list_0[ad] -= A_list[A_summary(N_phi, j1, j2, j3, j6, j5, j4)];
							A_list_0[ad] += A_list[A_summary(N_phi, j1, j3, j2, j6, j5, j4)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j2, j3, j1, j6, j5, j4)];
							A_list_0[ad] += A_list[A_summary(N_phi, j2, j1, j3, j6, j5, j4)];
							A_list_0[ad] -= A_list[A_summary(N_phi, j3, j1, j2, j6, j5, j4)];
							A_list_0[ad] += A_list[A_summary(N_phi, j3, j2, j1, j6, j5, j4)];
						}
					}
				}
			}
		}
	}
}

inline size_t A_summary_0(int N_phi, int j1, int j2, int j3, int j4, int j5, int j6)
{
	int j16, j43, j62, j24;
	size_t out = 0;
	size_t s = 2*N_phi - 1;

	j16 = (j1 - j6 + N_phi - 1);
	j43 = (j4 - j3 + N_phi - 1);
	j62 = (j6 - j2 + N_phi - 1);
	j24 = (j2 - j4 + N_phi - 1);

	out = s*s*s*j16 + s*s*j43 + s*j62 + j24;

	return out;
}

inline size_t A_summary(int N_phi, int j1, int j2, int j3, int j4, int j5, int j6)
{
	int j16, j43, j62, j24;
	size_t out = 0;
	size_t s = 2*N_phi - 1;

	j16 = (j1 - j6);
	j43 = (j4 - j3);
	j62 = (j6 - j2);
	j24 = (j2 - j4);
	//printf("%d,%d,%d,%d\t",j16+ N_phi - 1,j43+ N_phi - 1,j62+ N_phi - 1,j24+ N_phi - 1);
	j16 += ((j16<0)?s:(N_phi-1));
	j43 += ((j43<0)?s:(N_phi-1));
	j62 += ((j62<0)?s:(N_phi-1));
	j24 += ((j24<0)?s:(N_phi-1));
	//printf("%d,%d,%d,%d\n",j16,j43,j62,j24);

	out = s*s*s*j16 + s*s*j43 + s*j62 + j24;

	return out;
}
