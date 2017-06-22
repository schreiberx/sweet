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
scons --threading=omp --unit-test=test_rexi --gui=disable
EXEC="./build/test_rexi_*_release -N 64 --rexi-l=11 --rexi-normalization 0 --timestepping-method=ln_erk --timestepping-order 4"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
