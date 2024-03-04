

#pragma once

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <string.h>
#include <array>
#include <omp.h>
#include <math.h>
#include <complex>
#include <algorithm>
#include <fstream>


#define __DATA_DOUBLE
#ifdef __DATA_DOUBLE
//using MKL_Complex16 = double;
using Data_dtype = double;
#else
using MKL_Complex16 = std::complex<double>;
using Data_dtype = std::complex<double>;
#endif

using MKL_INT = unsigned long;
//#include<mkl.h>

using Data_basis = unsigned char;
//#include<mkl.h>
//#include<mkl_spblas.h>

//#define __PRINT_BASIS
#ifdef __PRINT_BASIS
#define __PRINT_BASIS_DETAIL
#define __PRINT_Z2INDEX
#endif

#define __TIMER
#ifdef __TIMER
#include "time/time.cc"
#endif

#ifndef __PRINTCG
//#define __PRINTCG
#endif

//--------------------------------------------------------------------------------------------------------------------
//------------------------- declaration ------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------

//#define bool	_Bool
//#define true	1
//#define false	0
#define ULIMAX 0xFFFFFFFFFFFFFFFF

#include <boost/math/special_functions/beta.hpp>
#include <boost/format.hpp>

#include "Basis/read.h"
#include "Basis/Basis.h"