
void cal2BodyAlist(struct CG2 *cg, double *pseudo, double *A_list2, int num_thread)
{
	double *cgt = NULL;
	long int u = (size_t)cg->u;
	double temp = 0.0;
	memset(A_list2, 0, u*u*u*u);

	# pragma omp parallel for num_threads(num_thread)
	for(long int m1=0;m1<u;++m1)
		for(long int m2=0;m2<u;++m2)
			for(long int m3=0;m3<u;++m3)
				for(long int m4=0;m4<u;++m4)
					for(long int j=0;j<u;++j)
					{
						temp = u/(2*j+1);
						temp *= temp;
						temp *= cg->cg[j][(u-1-m1)*u+u-1-m4]*cg->cg[j][(u-1-m3)*u+u-1-m2]*pseudo[j];
						temp *= cg->cg[j][0]*cg->cg[j][0];
						temp *= ((m1+m3)&1)?(-1.0):1.0;
						A_list2[m1*u*u*u + m2*u*u + m3*u + m4] += temp;
					}
	return;
}

