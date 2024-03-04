
MKL_Complex16 z_dot_c(const struct zvec* x, const struct zvec* y);

inline double z_norm(const struct zvec* x);

double z_normalization(struct zvec* x);

MKL_Complex16 z_orthogonlization(struct zvec* y, const struct zvec* x); // y = y - <x,y>x return <x,y>;

inline void z_axpy(struct zvec* y, const struct zvec* x, MKL_Complex16 a) ;

void z_copy(const struct zvec* x, struct zvec* y);

void print_z_vec(const struct zvec x);

void print_z(const MKL_Complex16 x);