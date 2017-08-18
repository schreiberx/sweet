#! /bin/bash


echo "***********************************************"
echo "Running tests for REXI with simple PDE"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST REXI PDE without halving"
echo "***********************************************"
#make clean
scons --threading=omp --unit-test=test_rexi_pde --gui=disable
EXEC="./build/test_rexi_pde_*_release -N 64 --rexi-l=11 --rexi-normalization 0 --rexi-half 1 --rexi-h=0.2 --rexi-m=512"
echo "$EXEC"
$EXEC || exit


echo
echo "***********************************************"
echo "TEST REXI PDE with halving"
echo "***********************************************"
#make clean
scons --threading=omp --unit-test=test_rexi_pde --gui=disable
EXEC="./build/test_rexi_pde_*_release -N 64 --rexi-l=11 --rexi-normalization 0 --rexi-half 0 --rexi-h=0.2 --rexi-m=512"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
