#! /bin/bash


echo "***********************************************"
echo "Running tests for REXI"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST ADVECTION: convergence in space (release) $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_rexi --gui=disable
EXEC="./build/test_rexi_libfft_omp_gnu_release"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
