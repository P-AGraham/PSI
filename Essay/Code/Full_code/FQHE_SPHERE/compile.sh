#!/bin/bash

# FQHE_SPHERE
cd ./core/Hilbertspace/SLFermion/

#icc cFQHE_Sphere.c -fPIC -fopenmp -shared -o fermionSL.so -O3

cd ../QHFMFermion/
icpc QHFM_fermion.c -fPIC -fopenmp -shared -o QHFM_fermion_real.so -O3
#icpc QHFM_fermion.c -fopenmp -o test -O3

cd ../../Lanczos/
icc My_spmv.c -fPIC -shared -o My_spmv.so -fopenmp
#icc Lanczos.c -fPIC -shared -o Lanczos.so  -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -DMKL_ILP64  -I"${MKLROOT}/include" -fopenmp
