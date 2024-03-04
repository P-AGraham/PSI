#!/bin/bash
FLAG='-L${MKLROOT}/lib/intel64/ -L${BOOST_LIB}/ -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl  -DMKL_ILP64  -I ${MKLROOT}/include/ -I ${BOOST_INCLUDE}/'

g++ -I ${BOOST_INCLUDE} QHFM_run.cc -fopenmp -o test -O2 ${FLAG} -std=c++17

if [ $# -gt 0 ];then
	echo "debug mode\n"
	g++ -I ${BOOST_INCLUDE} QHFM_run.cc -fopenmp -o test-g -O2 ${FLAG} -std=c++17 -g
fi
