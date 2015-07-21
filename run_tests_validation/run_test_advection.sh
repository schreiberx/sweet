#! /bin/bash


echo "***********************************************"
echo "Running tests for Spectral operations"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST ADVECTION (release) $X"
echo "***********************************************"
make clean
scons --unit-test=test_advection --gui=disable --spectral-space=disable --mode=release --spectral-dealiasing=disable
EXEC="./build/test_advection_gnu_release -X 1 -Y 1 -a 0.1 -N 32 -c 2 -d 1 -G 0 -C 0.2 -H 0 -R 1"
echo "$EXEC"
$EXEC || exit


echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
