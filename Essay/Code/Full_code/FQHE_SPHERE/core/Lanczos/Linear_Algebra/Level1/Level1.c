
MKL_Complex16 z_dot_c(const struct zvec* x, const struct zvec* y) // <x, y>
{
	MKL_Complex16 z;
	cblas_zdotc_sub(x->len, x->vec, 1, y->vec, 1, &z);
	return z;
}

inline double z_norm(const struct zvec* x)
{
	return cblas_dznrm2(x->len, x->vec, 1);
}

double z_normalization(struct zvec* x)
{
	double norm_z = z_norm(x);
	cblas_zdscal(x->len, 1.0/norm_z, x->vec, 1);
	return norm_z;
}

MKL_Complex16 z_orthogonlization(struct zvec* y, const struct zvec* x) 
// cal y = y-<x,y>*x and return <x,y>
{
	MKL_Complex16 xy = -z_dot_c(x, y);
	cblas_zaxpy(x->len, &xy, x->vec, 1, y->vec, 1);
	return -xy;
}

inline void z_axpy(struct zvec* y, const struct zvec* x, MKL_Complex16 a) 
// cal y = y-<x,y>*x and return <x,y>
{
	cblas_zaxpy(x->len, &a, x->vec, 1, y->vec, 1);
	return;
}

void z_copy(const struct zvec* x, struct zvec* y)
{
	cblas_zcopy(x->len, x->vec, 1, y->vec, 1);
	y->len = x->len;
	return;
}

void print_z(const MKL_Complex16 x)
{
	printf("%e+%eI",creal(x),cimag(x));
	return;
}

void print_z_vec(const struct zvec x)
{
	for(size_t i=0;i<x.len;++i)
	{
		print_z(x.vec[i]);printf("\t");
		if(i%5==4)
			printf("\n");
	}
}