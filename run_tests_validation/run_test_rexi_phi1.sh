#! /bin/bash


echo "***********************************************"
echo "Running tests for REXI"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST REXI without halving"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_rexi_phi1 --gui=disable
EXEC="./build/test_rexi_phi1_*_release -N 64 --rexi-l=11 --rexi-normalization 0 --rexi-half 1"
echo "$EXEC"
$EXEC || exit


echo
echo "***********************************************"
echo "TEST REXI with halving"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_rexi_phi1 --gui=disable
EXEC="./build/test_rexi_phi1_*_release -N 64 --rexi-l=11 --rexi-normalization 0 --rexi-half 0"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
