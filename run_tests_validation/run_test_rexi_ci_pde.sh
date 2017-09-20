#! /bin/bash


echo "***********************************************"
echo "Running tests for REXI with simple PDE"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST REXI CI PDE"
echo "***********************************************"
#make clean
scons --threading=omp --unit-test=test_rexi_ci_pde --gui=disable
EXEC="./build/test_rexi_ci_pde_*_release --rexi-method=ci --rexi-ci-primitive=circle --rexi-ci-mu=0 --rexi-ci-n=64 --rexi-ci-sx=25 --rexi-ci-sy=25 -v 5"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
