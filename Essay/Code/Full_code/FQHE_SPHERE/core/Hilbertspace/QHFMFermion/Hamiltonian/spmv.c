


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



void my_spmv_d(double* data, size_t* indices, size_t* indptr, double* vector, double* out, size_t dim, size_t n_threads)
{
    size_t i,j;

    # pragma omp parallel for num_threads(n_threads)
    for(i=0;i<dim;i++)
    {   
        out[i] *= 0;
        for(j=indptr[i];j<indptr[i+1];j++)
        {
            out[i] += data[j] * vector[indices[j]];
        }
    }
}