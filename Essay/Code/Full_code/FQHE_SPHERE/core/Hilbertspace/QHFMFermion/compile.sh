#!/bin/bash
FLAG='-L${MKLROOT}/lib/intel64/ -L${BOOST_LIB}/ -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl  -DMKL_ILP64  -I ${MKLROOT}/include/ -I ${BOOST_INCLUDE}/'

icpc -I ${BOOST_INCLUDE} RES.cc -fopenmp -o test -O3 ${FLAG} -std=c++17
