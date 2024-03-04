#include <iostream>
#include <unordered_map>
#include <string.h>
#include <array>
#include <omp.h>
#include <math.h>
#include <complex>
#define MKL_INT unsigned int
#define MKL_Complex16 complex<double>
#define Data_dtype complex<double>
using namespace std;

//#define I complex<double> {0,1};

//#include "Basis/Basis.h"
//#include "Hamiltonian/Hamiltonian.h"
#include "DLFermion.h"

//#include "Hamiltonian/Hamiltonian.cpp"
//#include "Basis/Basis.cpp"


int main()
{
	init(10,2,2,6);
	cout << A << endl;
	create_k1_space(0);
	create_k2_space(0);
	size_t k2_dim = get_k2_dim();
	cout << "k2_dim = " << k2_dim << endl;
	unsigned int *indptr = new unsigned int [k2_dim+1];
	cout << A.count_2symm(indptr) << endl;
	unsigned int* indices = new unsigned int [indptr[k2_dim]];
	complex<double>* data = new complex<double> [indptr[k2_dim]];
	complex<double>* A1_list = new complex<double> [(2*A.get_N_phi()-1) * (2*A.get_N_phi()-1)];
	complex<double>* A2_list = new complex<double> [(2*A.get_N_phi()-1) * (2*A.get_N_phi()-1)];
	complex<double>* A12_list = new complex<double> [(2*A.get_N_phi()-1) * (2*A.get_N_phi()-1)];
	cout << hamiltonian_generate(A1_list, A2_list, A12_list, indptr, indices, data) << endl;
	create_k1_space(0);
	create_k2_space(0);
	cout << "k2_dim = " << get_k2_dim() << endl;
	complex<double>* ket_k1 = new complex<double> [A.get_k1_dim()];
	complex<double>* out = new complex<double> [binomial(3*A.get_N_phi(), A.get_N_up()+A.get_N_down())];
	A.to_3layers_trans(ket_k1, out);
	delete [] indptr,indices,A1_list,A2_list,A12_list,ket_k1,out;
	return 0;
}

















































/*
	DLFermion *b = new DLFermion(30, 10, 0, 6);
	DLFermion & a = *b;
	cout << a << endl;
	a.create_k1_basis(0);
	a.create_k2_basis(0);
	size_t k2_dim = a.get_k2_dim();
	cout << "k2_dim = " << k2_dim << endl;
	unsigned int *indptr = new unsigned int [k2_dim+1];
	cout << a.count_2symm(indptr) << endl;
	unsigned int* indices = new unsigned int [indptr[k2_dim]];
	complex<double>* data = new complex<double> [indptr[k2_dim]];
	complex<double>* A1_list = new complex<double> [(2*a.get_N_phi()-1) * (2*a.get_N_phi()-1)];
	complex<double>* A2_list = new complex<double> [(2*a.get_N_phi()-1) * (2*a.get_N_phi()-1)];
	complex<double>* A12_list = new complex<double> [(2*a.get_N_phi()-1) * (2*a.get_N_phi()-1)];
	cout << a.generate_2symm(A1_list, A2_list, A12_list, indptr, indices, data) << endl;
	delete [] indptr,indices,A1_list,A2_list,A12_list;

	a.reset(16,4,4,6);
	cout << a << endl;
	a.create_k1_basis(0);
	a.create_k2_basis(0);
	k2_dim = a.get_k2_dim();
	cout << "k2_dim = " << k2_dim << endl;
	indptr = new unsigned int [k2_dim+1];
	cout << a.count_2symm(indptr) << endl;
	indices = new unsigned int [indptr[k2_dim]];
	data = new complex<double> [indptr[k2_dim]];
	A1_list = new complex<double> [(2*a.get_N_phi()-1) * (2*a.get_N_phi()-1)];
	A2_list = new complex<double> [(2*a.get_N_phi()-1) * (2*a.get_N_phi()-1)];
	A12_list = new complex<double> [(2*a.get_N_phi()-1) * (2*a.get_N_phi()-1)];
	cout << a.generate_2symm(A1_list, A2_list, A12_list, indptr, indices, data) << endl;
	delete b;
	delete [] indptr,indices,A1_list,A2_list,A12_list;
*/