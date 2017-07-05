#! /bin/bash


echo "***********************************************"
echo "Running tests for REXI"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST REXI: convergence in space (release) $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_rexi_ng --gui=disable
EXEC="./build/test_rexi_ng_*_release -N 64"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
